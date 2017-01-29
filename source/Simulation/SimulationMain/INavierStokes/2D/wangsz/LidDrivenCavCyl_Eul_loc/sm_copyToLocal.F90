! Test the array 
! Shizhao Wang, Jan 19, 2015

!#define USE_1D_ARRAY
#define USE_RESHAPE
 
subroutine sm_copyToLocal(gridVar,Array,nx,ny,nz,nv,flag)

  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getBlkIndexLimits, Grid_getListOfBlocks, Grid_getBlkCornerID

  implicit none
#include "constants.h"
#include "Flash.h"

  integer,intent(IN) :: gridVar,flag
  integer, intent(IN) :: nx, ny, nz, nv
#ifdef USE_1D_ARRAY
  real, dimension(nx*ny*nz) :: Array
#else
  real, dimension(nx,ny,nz),intent(OUT):: Array
#endif

  real, dimension(:,:,:,:), pointer :: solnData
  real, dimension(:,:,:,:), pointer :: facexData, faceyData, facezData

  integer :: blockID, i, j, k, n, iblk, is, fanout

#ifdef USE_RESHAPE
  real, allocatable, dimension(:,:,:,:) :: tmp
  integer :: newShape(4), order(4)
#endif

     !This was only applicable to UG, now also works for 1D PM - KW
        blockID = 1
#ifdef USE_1D_ARRAY
!        write(*,*) 'Copy to local, 1D array'
        if(flag == CENTER) then
          call Grid_getBlkPtr(blockID,solnData,CENTER)
          n = 0
          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                n = n + 1
                Array(n) = solnData(gridVar,i,j,k)
              enddo
            enddo
          enddo
          call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        elseif(flag == FACEX) then
          call Grid_getBlkPtr(blockID,facexData,FACEX)
          n = 0
          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                n = n + 1
                Array(n) = facexData(gridVar,i,j,k)
              enddo
            enddo
          enddo
          call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        elseif(flag == FACEY) then
          call Grid_getBlkPtr(blockID,faceyData,FACEY)
          n = 0
          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                n = n + 1
                Array(n) = faceyData(gridVar,i,j,k)
              enddo
            enddo
          enddo
          call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        elseif(flag == FACEZ) then
          call Grid_getBlkPtr(blockID,facezData,FACEZ)
          n = 0
          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                n = n + 1
                Array(n) = facezData(gridVar,i,j,k)
              enddo
            enddo
          enddo
          call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
        else
          write(*,*) 'Error in copy to local'
          stop
        endif
#else
#ifdef USE_RESHAPE
        allocate(tmp(nx,ny,nz,nv))
        newShape(1) = nx
        newShape(2) = ny
        newShape(3) = nz
        newShape(4) = nv
        order(1)    = 4
        order(2)    = 1
        order(3)    = 2
        order(4)    = 3

        if(flag == CENTER) then
          call Grid_getBlkPtr(blockID,solnData,CENTER)
          tmp = reshape(solnData,newShape, ORDER=order)
          Array(:,:,:) = tmp(:,:,:,gridVar)
          call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        elseif(flag == FACEX) then
          call Grid_getBlkPtr(blockID,facexData,FACEX)
          tmp = reshape(facexData,newShape, ORDER=order)
          Array(:,:,:) = tmp(:,:,:,gridVar)
          call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        elseif(flag == FACEY) then
          call Grid_getBlkPtr(blockID,faceyData,FACEY)
          tmp = reshape(faceyData,newShape, ORDER=order)
          Array(:,:,:) = tmp(:,:,:,gridVar)
          call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        elseif(flag == FACEZ) then
          call Grid_getBlkPtr(blockID,facezData,FACEZ)
          tmp = reshape(facezData,newShape, ORDER=order)
          Array(:,:,:) = tmp(:,:,:,gridVar)
          call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
        else
          write(*,*) 'Error in copy to local'
          stop
        endif
#else
        if(flag == CENTER) then
          call Grid_getBlkPtr(blockID,solnData,CENTER)
          Array(:,:,:) = solnData(gridVar,:,:,:)
          call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        elseif(flag == FACEX) then
          call Grid_getBlkPtr(blockID,facexData,FACEX)
          Array(:,:,:) = facexData(gridVar,:,:,:)
          call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        elseif(flag == FACEY) then
          call Grid_getBlkPtr(blockID,faceyData,FACEY)
          Array(:,:,:) = faceyData(gridVar,:,:,:)
          call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        elseif(flag == FACEZ) then
          call Grid_getBlkPtr(blockID,facezData,FACEZ)
          Array(:,:,:) = facezData(gridVar,:,:,:)
          call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
        else
          write(*,*) 'Error in copy to local'
          stop
        endif
#endif
#endif

  return
end subroutine sm_copyToLocal
