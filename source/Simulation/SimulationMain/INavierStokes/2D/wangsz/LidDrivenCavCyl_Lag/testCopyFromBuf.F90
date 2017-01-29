! Test the array 
! Shizhao Wang, Jan 19, 2015

subroutine testCopyFromBuf(gridVar,pfft_inArray)

  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getBlkIndexLimits, Grid_getListOfBlocks, Grid_getBlkCornerID

  implicit none
#include "constants.h"
#include "Flash.h"
  integer,intent(IN) :: gridVar
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC),intent(OUT):: pfft_inArray

  real, dimension(:,:,:,:), pointer :: solnData

  integer :: blockID, i, j, k, n, iblk, is, fanout


     !This was only applicable to UG, now also works for 1D PM - KW
        blockID = 1
        call Grid_getBlkPtr(blockID,solnData,CENTER)

        solnData(gridVar,:,:,:) = pfft_inArray(:,:,:)

        call Grid_releaseBlkPtr(blockID,solnData,CENTER)


end subroutine testCopyFromBuf
