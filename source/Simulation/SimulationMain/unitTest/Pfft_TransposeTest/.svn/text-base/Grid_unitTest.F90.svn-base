!!****if* source/Simulation/SimulationMain/unitTest/Pfft_TransposeTest/Grid_unitTest
!!
!! NAME
!!
!!  Grid_unitTest
!!
!! SYNOPSIS
!!
!!  Grid_unitTest(integer, intent(in):: fileUnit,
!!                logical, intent(inout)::perfect  )
!!
!! DESCRIPTION
!!
!!  This unit test checks that DENS_VAR contains the same as PDEN_VAR.
!!
!! ARGUMENTS
!!
!!  fileUnit - open f90 write unit
!!  perfect - returns a true if the test passed, false otherwise
!!
!! NOTES
!!
!!***

subroutine Grid_unitTest(fileUnit,perfect)
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getListOfBlocks
  implicit none
#include "Flash.h"
#include "constants.h"
  integer, intent(in)           :: fileUnit ! Output to file
  logical, intent(inout)        :: perfect  ! Flag to indicate errors

  real, pointer, dimension(:,:,:,:) :: solnData
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MAXBLOCKS) :: blockList
  integer :: blockCount, blockID, blk
  integer :: i, j, k

  call Grid_getListOfBlocks(LEAF, blockList, blockCount)
  do blk = 1, blockCount
     blockID = blockList(blk)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
           do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
              if ( int(solnData(DENS_VAR,i,j,k)) /= &
                   int(solnData(PDEN_VAR,i,j,k)) ) then
                 perfect = .false.
              end if
           end do
        end do
     end do
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  end do
end subroutine Grid_unitTest
