subroutine sim_compare_with_analytic_solution(mincellsize, globalsize, bbxlims_global, ismooth)

  use Simulation_data
  use Grid_interface, ONLY: Grid_getListOfBlocks,Grid_getBlkIndexLimits, &
                            Grid_getCellCoords, Grid_getPointData, Grid_putPointData, &
                            Grid_getDomainBoundBox
  use Grid_data, ONLY: gr_meshComm
  use IO_interface, ONLY: IO_writePlotfile
  use tree, ONLY : lrefine_max
  use Logfile_interface, ONLY: Logfile_stampMessage

  implicit none
#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"

  real, intent(IN) :: mincellsize
  integer, dimension(MDIM), intent(IN) :: globalSize
  real, dimension(2), intent(IN) :: bbxlims_global
  integer, intent(IN) :: ismooth

  integer :: blkCount, blockID, iblk, ierr
  integer, dimension(MAXBLOCKS) :: blkList
  integer :: i, j, k, d
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, allocatable, dimension(:) :: iCoords, jCoords, kCoords
  integer, dimension(MDIM) :: cell
  integer :: isize, jsize, ksize
  real, dimension(2, MDIM) :: bbox
  real, dimension(NDIM) :: domain_ctr, gwidth
  real, dimension(3) :: mpiar_buf_send, mpiar_buf_recv
  real :: flsm_predicted, flsm_actual, flam_actual, auxv, arg1, arg2, arg3, norm, invrt2pi
  real :: sum_orig, sum_convolved, sum_predicted
  character(len=128) :: msgbuf

  ! Need the domain center to compute predicted model
  call Grid_getDomainBoundBox(bbox)
  invrt2pi = 1.0 / sqrt(2.0 * PI)  ! 1/sqrt(2*pi)
  norm = 1.0
  do d = 1, NDIM

    domain_ctr(d) = 0.5 * ( bbox(2, d) + bbox(1,d) )
    gwidth(d) = gwidth_rel(d) * 0.5 * ( bbox(2, d) - bbox(1,d) )
    gwidth(d) = sqrt(gwidth(d)**2 + (mincellsize*sim_smooth_radius)**2) ! We convolved it, remember?
    norm = norm * invrt2pi / gwidth(d)

  end do

  ! Get list of blocks at maximum refinement level
  call Grid_getListOfBlocks(REFINEMENT,blkList, blkCount, lrefine_max)

  ! Figure out storage
  blockID = blklist(1)
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  isize = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
  jsize = blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
  ksize = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1
  allocate(iCoords(isize))
  allocate(jCoords(jsize))
  allocate(kCoords(ksize))

  sum_orig = 0.0
  sum_convolved = 0.0
  sum_predicted = 0.0
! Loop over blocks, do the comparison thing
  do iblk = 1, blkCount
  
    blockID = blklist(iblk)
    call Grid_getCellCoords(IAXIS,blockID,CENTER,.false.,iCoords,isize)
    call Grid_getCellCoords(JAXIS,blockID,CENTER,.false.,jCoords,jsize)
    call Grid_getCellCoords(KAXIS,blockID,CENTER,.false.,kCoords,ksize)
    
    do i = 1, isize
      cell(1) = i
      arg1 = ( ( iCoords(i) - domain_ctr(1) ) / gwidth(1) )**2
      
      do j = 1, jsize
        cell(2) = j
        arg2 = ( ( jCoords(j) - domain_ctr(2) ) / gwidth(2) )**2
        
        do k = 1, ksize
          cell(3) = k
          arg3 = ( ( kCoords(k) - domain_ctr(3) ) / gwidth(3) )**2
	  
          call Grid_getPointData(blockId, CENTER, FLSM_MSCALAR, INTERIOR, cell, flsm_actual)
          call Grid_getPointData(blockId, CENTER, FLAM_MSCALAR, INTERIOR, cell, flam_actual)

          flsm_predicted = norm * exp(-0.5 * (arg1 + arg2 + arg3) )
          auxv = abs( flsm_actual - flsm_predicted )
          
          call Grid_putPointData(blockId, CENTER, AUXV_MSCALAR, INTERIOR, cell, auxv)
          call Grid_putPointData(blockId, CENTER, DENS_MSCALAR, INTERIOR, cell, flsm_predicted)
		  if ( iCoords(i) .ge. bbxlims_global(1) .and. iCoords(i) .le. bbxlims_global(2) ) then 
            sum_orig = sum_orig + flam_actual
            sum_convolved = sum_convolved + flsm_actual
            sum_predicted = sum_predicted + flsm_predicted
          endif
	  
        end do
      end do
    end do

  end do

  mpiar_buf_send(1) = sum_orig
  mpiar_buf_send(2) = sum_convolved
  mpiar_buf_send(3) = sum_predicted

  call MPI_AllReduce(mpiar_buf_send, mpiar_buf_recv, 4, FLASH_REAL, MPI_SUM, gr_meshComm, ierr)

  sum_orig = mpiar_buf_recv(1) * mincellsize**3
  sum_convolved = mpiar_buf_recv(2) * mincellsize**3
  sum_predicted = mpiar_buf_recv(3) * mincellsize**3

  write(msgbuf, '("Smoothing length ", i2, ":")') ismooth + 1
  call Logfile_stampMessage(msgbuf)
  write(msgbuf, '("Normalization of original distribution = ", 1pE12.5)') sum_orig
  call Logfile_stampMessage(msgbuf)
  write(msgbuf, '("Normalization of convolved distribution = ", 1pE12.5)') sum_convolved
  call Logfile_stampMessage(msgbuf)
  write(msgbuf, '("Normalization of predicted distribution = ", 1pE12.5)') sum_predicted
  call Logfile_stampMessage(msgbuf)

  deallocate(iCoords)
  deallocate(jCoords)
  deallocate(kCoords)

end subroutine sim_compare_with_analytic_solution
