! Test the array 
! Shizhao Wang, Jan 19, 2015

subroutine sm_copyToGlobal(gridVar,Array,nx,ny,nz,flag)

use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_getBlkIndexLimits, Grid_getListOfBlocks, Grid_getBlkCornerID
 
implicit none
#include "constants.h"
#include "Flash.h"

integer,intent(IN) :: gridVar,flag
integer, intent(IN) :: nx, ny, nz
real, dimension(nx,ny,nz),intent(IN):: Array
real, dimension(:,:,:,:), pointer :: solnData
real, dimension(:,:,:,:), pointer :: facexData, faceyData, facezData
integer :: blockID, i, j, k, n, iblk, is, fanout

!This was only applicable to UG, now also works for 1D PM - KW
  blockID = 1
    
  if(flag == CENTER) then
    call Grid_getBlkPtr(blockID,solnData,CENTER)
    solnData(gridVar,:,:,:) = Array(:,:,:)
    call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  elseif(flag == FACEX) then
    call Grid_getBlkPtr(blockID,facexData,FACEX)
    facexData(gridVar,:,:,:) = Array(:,:,:)
    call Grid_releaseBlkPtr(blockID,facexData,FACEX)
  elseif(flag == FACEY) then
    call Grid_getBlkPtr(blockID,faceyData,FACEY)
    faceyData(gridVar,:,:,:) = Array(:,:,:)
    call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
  elseif(flag == FACEZ) then
    call Grid_getBlkPtr(blockID,facezData,FACEZ)
    facezData(gridVar,:,:,:) = Array(:,:,:)
    call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
  else
    write(*,*) 'Error in copy to local'
    stop
  endif
  
  return
  end subroutine sm_copyToGlobal
