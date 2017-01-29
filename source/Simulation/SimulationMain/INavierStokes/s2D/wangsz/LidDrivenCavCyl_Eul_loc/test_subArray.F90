! Test the time used in call the subroutine with different type of
! arraies
! Shizhao Wang
! Jan 11, 2015

#include "constants.h"
#include "Flash.h"
#include "ImBound.h"

        subroutine test_subArray(a,b)
        !subroutine test_subArray(a,b,nx,ny,nz)
        use Timers_interface, ONLY : Timers_start, Timers_stop
        implicit none
        real :: a(GRID_IHI_GC+1,GRID_JHI_GC  ,GRID_KHI_GC)
        real :: b(GRID_IHI_GC,GRID_JHI_GC+1  ,GRID_KHI_GC)
        !integer :: nx, ny, nz
        !real :: a(nx+1, ny, nz)
        !real :: b(nx, ny+1, nz)
        !real :: b(NFACE_VARS, GRID_IHI_GC,GRID_JHI_GC+1  ,GRID_KHI_GC)
        !real :: b(NFACE_VARS, GRID_IHI_GC,GRID_JHI_GC+1  ,GRID_KHI_GC)
        integer :: i, j, k

        call Timers_start('onebyone')
        do k = 1, GRID_KHI_GC
          do j = 1, GRID_JHI_GC
            do i = 1, GRID_IHI_GC+1
              a(i,j,k) = 1.0
            enddo
          enddo
        enddo

        do k = 1, GRID_KHI_GC
          do j = 1, GRID_JHI_GC+1
            do i = 1, GRID_IHI_GC
              b(i,j,k) = 2.0
            enddo
          enddo
        enddo
        call Timers_stop('onebyone')

        call Timers_start('oneCol')
        a = 1.0
        b = 2.0
        call Timers_stop('oneCol')

        call Timers_start('oneCol2')
        a(:,:,:) = 1.0
        b(:,:,:) = 2.0
        call Timers_stop('oneCol2')

        return
        endsubroutine test_subArray

