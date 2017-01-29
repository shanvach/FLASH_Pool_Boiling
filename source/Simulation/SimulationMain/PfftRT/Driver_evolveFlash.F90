subroutine Driver_evolveFlash()

  use Driver_data, ONLY: dr_globalMe, dr_globalNumProcs
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getSimTime
  use Grid_data, ONLY: gr_meshComm
  use Grid_interface, ONLY: Grid_getListOfBlocks, Grid_getBlkBoundBox, &
                            Grid_getGlobalIndexLimits, Grid_getDomainBoundBox, &
                            Grid_getMinCellSize, &
                            Grid_pfftInit, Grid_pfft, Grid_pfftGetIndexLimits, &
                            Grid_pfftMapToInput, Grid_pfftMapFromOutput, &
                            Grid_pfftFinalize, Grid_getBlkPtr, Grid_releaseBlkPtr, &
                            Grid_markRefineSpecialized, Grid_updateRefinement, &
                            Grid_getBlkIndexLimits, Grid_putPointData
  use gr_pfftData, ONLY : pfft_inLen, pfft_outLen, pfft_usableProc
  use IO_interface, ONLY: IO_writePlotfile
  use Simulation_data
!  use Simulation_interface, ONLY : sim_convolve, sim_compare_with_analytic_solution
  use Simulation_interface, ONLY : sim_convolve
  use tree, ONLY : lrefine_max, surr_blks
  use Logfile_interface, ONLY: Logfile_stampMessage

  implicit none
#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"
#include "Pfft.h"

  interface
	subroutine make_convolution_gaussian(region_bb, sigma, mincellsize)
		implicit none
		real, intent(IN), dimension(2, MDIM) :: region_bb
		real, intent(IN) :: sigma, mincellsize
	end subroutine
  end interface

  real, dimension(2, MDIM) :: bb, region_bb
  real, dimension(2) :: bbxlims, bbxlims_global
  real :: mincellsize, sigma
  real, allocatable, dimension(:) :: data_dir, data_trans, data_conv, kernel_dir, kernel_trans
  real, allocatable, dimension(:) :: data_trans_buf
  real, pointer, dimension(:,:,:,:) :: solnData

  integer :: i, ierr, size_dir, size_trans, isfirst, ismooth, var_idx
  integer :: blkCount
  integer, dimension(MAXBLOCKS) :: blkList, blkListForCopy
  integer, dimension(MDIM) :: globalSize, localsize, ttype
  integer,dimension(LOW:HIGH,MDIM) :: limits_dir, limits_trans
  real :: simulation_time
  integer :: n_step

  integer :: iSizeGC, jSizeGC, kSizeGC
  integer :: iSize, jSize, kSize
  integer,dimension(LOW:HIGH,MDIM)::blkLimits,blkLimitsGC
  integer :: iblk, j, k, istat
  real :: flsm_val
  integer, dimension(MDIM) :: cell_position
  real, allocatable, dimension(:) :: xCenter
  real :: domainWidth
  real :: smooth_radius

  logical :: needMap, gridChanged

  character(len=128) :: msgbuf


  !############################################################################
  ! This routine is used to analyze smoothings of the flame surface in 
  ! the highest refinement region of the RT flame simulation. The region
  ! of highest refinement is assumed to be a rectangular prism that spans
  ! the y and z directions and occupies a fraction of the domain in the
  ! x direction. This assumption is valid for RT flame simulations that
  ! force this refinement pattern.
  !############################################################################

  !********************
  ! Get simulation time and write out initial data from original file
  !********************

  call Driver_getSimTime(simulation_time)
  call Driver_getNStep(n_step)

  isfirst=1
  
  call IO_writeIntegralQuantities(isfirst, simulation_time)

  isfirst=0

  !********************
  ! Copy FLAM to FLSM to fill in data outside of
  ! maximally refined region.
  !********************

  call Grid_getListOfBlocks(ALL_BLKS, blkListForCopy, blkCount)

  do i = 1, blkCount

     call Grid_getBlkPtr(blkListForCopy(i), solnData)

     solnData(FLSM_MSCALAR,:,:,:) = solnData(FLAM_MSCALAR,:,:,:)

     call Grid_releaseBlkPtr(blkListForCopy(i), solnData)

  end do

  !********************
  ! Figure out the original bounding box.
  !********************

  ! Establish the full region bounding box
  call Grid_getDomainBoundBox(region_bb)

  ! Get list of blocks at maximum refinement level
  call Grid_getListOfBlocks(REFINEMENT,blkList, blkCount, lrefine_max)

  ! Loop over the list of blocks and get all the bounding boxes.  Record
  ! the highest value of the left x bound and the lowest value of the
  ! right x bound.
  bbxlims(1) = region_bb(1, IAXIS)
  bbxlims(2) = region_bb(2, IAXIS)
  do i = 1, blkCount

    call Grid_getBlkBoundBox(blkList(i), bb)

    ! if the left edge of this block is a refinement boundary, and if its
    ! left boundary is the largest one so far, record it.
    if ( surr_blks(PROCNO, LEFT_EDGE, CENTER, CENTER, blklist(i)) == NONEXISTENT .and. &
         bb(1, IAXIS) > bbxlims(1) )         bbxlims(1) = bb(1, IAXIS)
    
    ! if the right edge of this block is a refinement boundary, and if its
    ! right boundary is the least one so far, record it.
    if ( surr_blks(PROCNO, RIGHT_EDGE, CENTER, CENTER, blklist(i)) == NONEXISTENT .and. &
         bb(2, IAXIS) < bbxlims(2) )         bbxlims(2) = bb(2, IAXIS)
  
  end do

  ! Get the global versions of bbxlims
  call MPI_AllReduce(bbxlims(1), bbxlims_global(1), 1, FLASH_REAL, MPI_MAX, gr_meshComm, ierr)
  call MPI_AllReduce(bbxlims(2), bbxlims_global(2), 1, FLASH_REAL, MPI_MIN, gr_meshComm, ierr)
  if ( dr_globalME == MASTER_PE .and. bbxlims_global(1) >= bbxlims_global(2) ) then
    call Driver_abortFLash("[Driver_evolveFlash]: No fully-refined slab found, bummer.")
  endif
  
  ! Now establish the region bounding box restricted to this fully-refined slab
  region_bb(1, IAXIS) = bbxlims_global(1)
  region_bb(2, IAXIS) = bbxlims_global(2)

  original_region_bb = region_bb

  !********************  

  !********************
  ! Now that the fully refined slab is found, 
  ! refine a buffer region of length bufFact*(ymax-ymin)
  ! above and below the fully refined region (i.e. make it bigger)
  !********************

  gridChanged = .true.

  do while (gridChanged) 

     call Grid_updateRefinement(1, simulation_time, gridChanged)

  end do 

  !********************

  !********************
  ! Relocate slab and reset values
  !********************

  ! Establish the full region bounding box
  call Grid_getDomainBoundBox(region_bb)


  ! Get list of blocks at maximum refinement level
  call Grid_getListOfBlocks(REFINEMENT,blkList, blkCount, lrefine_max)

  ! Loop over the list of blocks and get all the bounding boxes.  Record
  ! the highest value of the left x bound and the lowest value of the
  ! right x bound.
  bbxlims(1) = region_bb(1, IAXIS)
  bbxlims(2) = region_bb(2, IAXIS)
  do i = 1, blkCount

    call Grid_getBlkBoundBox(blkList(i), bb)

    ! if the left edge of this block is a refinement boundary, and if its
    ! left boundary is the largest one so far, record it.
    if ( surr_blks(PROCNO, LEFT_EDGE, CENTER, CENTER, blklist(i)) == NONEXISTENT .and. &
         bb(1, IAXIS) > bbxlims(1) )         bbxlims(1) = bb(1, IAXIS)

    ! if the right edge of this block is a refinement boundary, and if its
    ! right boundary is the least one so far, record it.
    if ( surr_blks(PROCNO, RIGHT_EDGE, CENTER, CENTER, blklist(i)) == NONEXISTENT .and. &
         bb(2, IAXIS) < bbxlims(2) )         bbxlims(2) = bb(2, IAXIS)

  end do

  ! Get the global versions of bbxlims
  call MPI_AllReduce(bbxlims(1), bbxlims_global(1), 1, FLASH_REAL, MPI_MAX,gr_meshComm, ierr)
  call MPI_AllReduce(bbxlims(2), bbxlims_global(2), 1, FLASH_REAL, MPI_MIN,gr_meshComm, ierr)
  if ( dr_globalME == MASTER_PE .and. bbxlims_global(1) >= bbxlims_global(2) ) then
    call Driver_abortFLash("[Driver_evolveFlash]: No fully-refined slab found,bummer.")
  endif

  ! Now establish the region bounding box restricted to this fully-refined slab
  region_bb(1, IAXIS) = bbxlims_global(1)
  region_bb(2, IAXIS) = bbxlims_global(2)

  !********************


  ! Get number of cells per dimension
  call Grid_getGlobalIndexLimits(globalSize)
  call Grid_getMinCellSize(mincellsize)
  globalsize(1) = ( bbxlims_global(2) - bbxlims_global(1) + 0.01 * mincellsize ) &
                   / mincellsize
  !********************
  ! Make the numerical Gaussian to be convolved with the data
  !********************

  !********************
  ! Make a little logfile noise
  !********************
  write(msgbuf, '("Slab Bounding Box = (",1pE12.5,",",1pE12.5,") ; Ncells = ",I9)') &
        bbxlims_global(1), bbxlims_global(2), globalsize(1)
  call Logfile_stampMessage(msgbuf)

  !********************
  ! Initialize pfft data structures.
  !********************

   needMap = ( NDIM > 1 )
   ttype(1) = PFFT_REAL2C_EXTEND
   ttype(2) = PFFT_COMPLEX
   ttype(3) = PFFT_COMPLEX
   call Grid_pfftInit(NDIM, needMap, globalsize, localsize, refinementLevel=NONEXISTENT, region_bndBox=region_bb, &
                      transformType = ttype )
 
  !********************
  ! Perform the convolution
  !********************
  if ( pfft_usableProc ) then
  
    size_dir = pfft_inLen(IAXIS) * pfft_inLen(JAXIS) * pfft_inLen(KAXIS)
    size_trans = 2 * pfft_outLen(IAXIS) * pfft_outLen(JAXIS) * pfft_outLen(KAXIS)
    allocate( data_dir(size_dir) )
    allocate( data_conv(size_dir) )
    allocate( kernel_dir(size_dir) )
    allocate( data_trans(size_trans) )
    allocate( data_trans_buf(size_trans))
    allocate( kernel_trans(size_trans) )
    
    call Grid_pfftGetIndexLimits(limits_dir,limits_trans)

    call Grid_pfftMapToInput(FLAM_MSCALAR,data_dir)
    call Grid_pfft(PFFT_FORWARD,data_dir,data_trans_buf)

  endif
 
  domainWidth = region_bb(2,JAXIS) - region_bb(1,JAXIS)

  do ismooth = sim_smooth_step_min, sim_smooth_step_max
		var_idx = FLSM_MSCALAR
		smooth_radius = 2.0**(real(ismooth)*sim_smooth_step_delta)
		sigma = mincellsize * smooth_radius
		write(msgbuf, '("Smoothing Length = ",1pE12.5)') sigma
  		call Logfile_stampMessage(msgbuf)
		call make_convolution_gaussian(region_bb, sigma, mincellsize)
		if ( pfft_usableProc ) then

			data_trans = data_trans_buf

			call Grid_pfftMapToInput(GAUS_MSCALAR,kernel_dir)
        
			call Grid_pfft(PFFT_FORWARD,kernel_dir,kernel_trans)
        
			call sim_convolve(data_trans, kernel_trans, limits_trans, globalsize)
        
			call Grid_pfft(PFFT_INVERSE, data_trans, data_conv)
			call Grid_pfftMapFromOutput(var_idx, data_conv)
		endif
!		call sim_compare_with_analytic_solution(mincellsize, globalsize, bbxlims_global, ismooth)

                !************************************
                ! Fix artifacts from convolution
                !************************************

                call Grid_getListOfBlocks(LEAF, blkList, blkCount)

                do iblk = 1, blkCount
             
                        call Grid_getBlkIndexLimits(blkList(iblk),blkLimits,blkLimitsGC)
                        iSizeGC = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
                        allocate(xCenter(iSizeGC),STAT=istat)
                        if (istat /= 0) call Driver_abortFlash("Cannot allocate xCenter in Driver_evolveFlash")

                        xCenter(:)=0.0
                        call Grid_getCellCoords(IAXIS,blkList(iblk), CENTER, .true.,xCenter,iSizeGC)

                        flsm_val = -1.0

                        do k = blkLimits(LOW,KAXIS) ,blkLimits(HIGH,KAXIS)
                           do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
                              do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
                                cell_position(1)=i
                                cell_position(2)=j
                                cell_position(3)=k

                                if (xCenter(i) < original_region_bb(1, IAXIS)) then
                                   flsm_val = 1.0
                                   call Grid_putPointData(blkList(iblk), &
                                        CENTER, FLSM_MSCALAR, EXTERIOR, & 
                                        cell_position, flsm_val) 
                               
                                else if (xCenter(i) > original_region_bb(2,IAXIS)) then
                                   flsm_val =  0.0
                                   call Grid_putPointData(blkList(iblk), &
                                        CENTER,FLSM_MSCALAR, EXTERIOR, &
                                        cell_position, flsm_val)
                                end if
                              end do
                           end do
                        end do

                        deallocate(xCenter, STAT=istat)
                        if (istat /= 0) call Driver_abortFlash("Cannot deallocate xCenter in Driver_evolveFlash")
                end do
  

		call IO_writeIntegralQuantities(isfirst, sigma)

                if (sim_writeflsmdata) then
		   call Grid_restrictAllLevels()
		   call IO_writePlotfile(.false.)
                endif
	
    end do

  if ( pfft_usableProc ) then
  
    deallocate(data_dir)
    deallocate(data_conv)
    deallocate(data_trans)
    deallocate(data_trans_buf)
    deallocate(kernel_dir)
    deallocate(kernel_trans)

  endif

  call Grid_pfftFinalize()

end subroutine

subroutine make_convolution_gaussian(region_bb, sigma, mincellsize)

	use Grid_interface, ONLY: Grid_getListOfBlocks,Grid_getBlkIndexLimits, &
							  Grid_getCellCoords, Grid_putPointData
	use tree, ONLY : lrefine_max
	use Simulation_data

	implicit none
#include "Flash.h"
#include "constants.h"

	real, intent(IN), dimension(2, MDIM) :: region_bb
	real, intent(IN) :: sigma, mincellsize

	real :: arg1, arg2, arg3, norm, invrt2pi, gaus, var
	integer :: blkCount, blockID, iblk
	integer, dimension(MAXBLOCKS) :: blkList
	integer :: i, j, k, d
	integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
	real, allocatable, dimension(:) :: iCoords, jCoords, kCoords
	integer :: isize, jsize, ksize
    integer, dimension(MDIM) :: cell

!	norm = (2.0 * PI)**(-1.5) / sigma**3
!	norm = (PI)**(-1.5) / sigma**3
        norm = ((PI)**(0.5))*sigma*erf((region_bb(2,JAXIS) - region_bb(1,JAXIS))*0.5/sigma)
        norm = norm*((PI)**(0.5))*sigma*erf((region_bb(2,KAXIS) - region_bb(1,KAXIS))*0.5/sigma)
        norm = norm*((PI)**(0.5))*sigma*erf((region_bb(2,IAXIS) - region_bb(1,IAXIS))*0.5/sigma)

	var = sigma * sigma
	
	domain_volume = 1.0
	do d = 1, MDIM
		gauss_ctr(d) = 0.5 * ( region_bb(2, d) + region_bb(1,d) )
		gauss_ctr_idx(d) = ( gauss_ctr(d) - region_bb(1,d) )/ mincellsize
		domain_volume = domain_volume * ( region_bb(2, d) - region_bb(1, d) )
	end do

	! Get list of blocks at maximum refinement level
	call Grid_getListOfBlocks(REFINEMENT, blkList, blkCount, lrefine_max)
	
	! Figure out storage
	blockID = blklist(1)
	call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
	isize = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
	jsize = blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
	ksize = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1
	allocate(iCoords(isize))
	allocate(jCoords(jsize))
	allocate(kCoords(ksize))

	do iblk = 1, blkCount
  
		blockID = blklist(iblk)
		call Grid_getCellCoords(IAXIS,blockID,CENTER,.false.,iCoords,isize)
		call Grid_getCellCoords(JAXIS,blockID,CENTER,.false.,jCoords,jsize)
		call Grid_getCellCoords(KAXIS,blockID,CENTER,.false.,kCoords,ksize)
		
		do i = 1, isize
			cell(1) = i
			arg1 = ( ( iCoords(i) - gauss_ctr(1) ) )**2
		
			do j = 1, jsize
			cell(2) = j
			arg2 = ( ( jCoords(j) - gauss_ctr(2) )  )**2
			
				do k = 1, ksize
					cell(3) = k
					arg3 = ( ( kCoords(k) - gauss_ctr(3) )  )**2
	
					!gaus = norm * exp(-0.5 * (arg1+arg2+arg3) / var)
					!gaus = norm * exp(-1.0 * (arg1+arg2+arg3) / var)
					gaus = exp(-1.0 * (arg1+arg2+arg3) / var)/norm
					
					call Grid_putPointData(blockId, CENTER, GAUS_MSCALAR, INTERIOR, cell, gaus)
				end do			
			end do
		end do
		
	end do

	deallocate(iCoords)
    deallocate(jCoords)
    deallocate(kCoords)

end subroutine
