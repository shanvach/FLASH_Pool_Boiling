!!****if* source/Grid/GridMain/gr_findMean
!!
!! NAME
!!  gr_findMean
!! 
!! SYNOPSIS
!!  gr_findMean(integer, intent(in)  :: iSrc,
!!                     integer, intent(in)  :: iType,
!!                     logical, intent(in)  :: bGuardcell
!!                        real, intent(out) :: mean)
!!
!! DESCRIPTION
!!  Calculates the mean of a function
!!
!!
!! ARGUMENTS
!!  iSrc -- the index (e.g. DENS_VAR) into the unk routine to calculate over
!!  iType -- the type of mean.  Valid values will be  (feel free to implement more)
!!     2 = arithmetic average
!!  bGuardcell -- logical indicating whether guard cells should be included in the calculation
!!  mean -- the requested mean
!!
!!
!!***
subroutine gr_findMean(iSrc, iType, bGuardcell, mean)
  
  use Driver_interface, ONLY: Driver_abortFlash
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
       Grid_getBlkPhysicalSize, Grid_getBlkIndexLimits, &
       Grid_getBlkPtr, Grid_releaseBlkPtr
  use Grid_data, ONLY : gr_meshComm
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer, intent(in) :: iSrc, iType
  logical, intent(in) :: bGuardcell
  real, intent(out) :: mean

  real :: localVolume, localSum, blockSum, sum

  real :: blockVolume, volume
  integer :: blkCount, lb, blockID, i, j, k, ierr
  integer :: ili, iui, jli, jui, kli, kui
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MAXBLOCKS) :: blkList
  real, dimension(:,:,:,:), pointer :: solnData
  real, dimension(GRID_IHI_GC) :: dx
  real, dimension(GRID_JHI_GC) :: dy
!!==============================================================================

  mean = 0.0
  localVolume = 0.
  localSum = 0.

  call Grid_getListOfBlocks(LEAF,blkList,blkCount)

  do lb = 1, blkCount

     blockID = blkList(lb)

     call Grid_getBlkPtr(blockID, solnData)
     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)

     blockSum = 0.
     blockVolume = 0.

     ili = blkLimits(LOW,IAXIS)
     iui = blkLImits(HIGH,IAXIS)
     jli = blkLimits(LOW,JAXIS)
     jui = blkLImits(HIGH,JAXIS)

     call Grid_getCellMetrics(IAXIS, lb, CENTER, .true., dx, GRID_IHI_GC)
     call Grid_getCellMetrics(JAXIS, lb, CENTER, .true., dy, GRID_JHI_GC) 

     do j = jli, jui
       do i = ili, iui
         blockSum = blockSum + solnData(iSrc,i,j,1) / ( dx(i) * dy(j) )
         blockVolume = blockVolume + 1 / (dx(i) * dy(j))
       enddo
     enddo

     localSum = localSum + blockSum 
     localVolume = localVolume + blockVolume

     call Grid_releaseBlkPtr(blockID, solnData)

  enddo

  call mpi_allreduce ( localSum, sum, 1, FLASH_REAL, MPI_SUM, gr_meshComm, ierr )
  call mpi_allreduce ( localVolume, volume, 1, FLASH_REAL, MPI_SUM, gr_meshComm, ierr )

  mean =  sum / volume

  return

end subroutine gr_findMean

