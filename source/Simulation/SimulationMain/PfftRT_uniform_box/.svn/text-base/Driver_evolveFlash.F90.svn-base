subroutine Driver_evolveFlash()

  use Driver_data, ONLY: dr_globalMe, dr_globalNumProcs
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY: gr_meshComm
  use Grid_interface, ONLY: Grid_getListOfBlocks, Grid_getBlkBoundBox, &
                            Grid_getGlobalIndexLimits, Grid_getDomainBoundBox, &
                            Grid_getMinCellSize, &
                            Grid_pfftInit, Grid_pfft, Grid_pfftGetIndexLimits, &
                            Grid_pfftMapToInput, Grid_pfftMapFromOutput, &
                            Grid_pfftFinalize
  use gr_pfftData, ONLY : pfft_inLen, pfft_outLen, pfft_usableProc
  use IO_interface, ONLY: IO_writeCheckpoint
  use Simulation_data
  use Simulation_interface, ONLY : sim_convolve, sim_compare_with_analytic_solution
  use tree, ONLY : lrefine_max, surr_blks
  use Logfile_interface, ONLY: Logfile_stampMessage

  implicit none
#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"
#include "Pfft.h"


  real, dimension(2, MDIM) :: bb, region_bb
  real, dimension(2) :: bbxlims, bbxlims_global
  real :: mincellsize, cwidth
  real, allocatable, dimension(:) :: data_dir, data_trans

  integer :: i, ierr, size_dir, size_trans
  integer :: blkCount
  integer, dimension(MAXBLOCKS) :: blkList
  integer, dimension(MDIM) :: globalSize, localsize, ttype
  integer,dimension(LOW:HIGH,MDIM) :: limits_dir, limits_trans

  logical :: needMap

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
  ! Figure out the bounding box.
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

  ! Get number of cells per dimension
  call Grid_getGlobalIndexLimits(globalSize)
  call Grid_getMinCellSize(mincellsize)
  globalsize(1) = ( bbxlims_global(2) - bbxlims_global(1) + 0.01 * mincellsize ) &
                   / mincellsize
  cwidth = box_hwidth * mincellsize

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
    allocate( data_trans(size_trans) )
    
    call Grid_pfftGetIndexLimits(limits_dir,limits_trans)
    call Grid_pfftMapToInput(FLAM_VAR,data_dir)
    call Grid_pfft(PFFT_FORWARD,data_dir,data_trans)
    call sim_convolve(data_trans, limits_trans, box_hwidth, globalsize)
    call Grid_pfft(PFFT_INVERSE, data_trans, data_dir)
    call Grid_pfftMapFromOutput(FLSM_VAR, data_dir)
    deallocate(data_dir)
    deallocate(data_trans)

  endif

  call Grid_pfftFinalize()

  call sim_compare_with_analytic_solution(mincellsize, globalsize, bbxlims_global)
  
  call IO_writeCheckpoint()
  

end subroutine
