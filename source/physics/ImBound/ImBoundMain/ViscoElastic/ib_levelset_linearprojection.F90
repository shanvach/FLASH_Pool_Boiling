!===============================================================================
!!
!! subroutine ib_levelset_linearprojection 
!!
!! extrapolating level set gradient in normal direction using 1st order upwind 
!!
!! args inout:: s       - levelset funtion to be projected
!! args in   :: so      - levelset at time n-1
!!           :: sn      - directional derivate of levelset 
!!           :: u,v     - normal vectors nx,ny
!!           :: dx,dy   - spacing in x and y directions
!!           :: ix1,ix2 - low,high block x indeces
!!           :: jy1,jy2 - low,high block y indeces 
!===============================================================================
#include "constants.h"
#include "Flash.h"

        subroutine ib_levelset_linearprojection(lmda,s,sn,u,v,ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz)
        implicit none
        !include 'mpif.h'
        real, dimension(:,:,:), intent(inout) :: s
        real, dimension(:,:,:), intent(in)    :: u,v,lmda
        real, dimension(:,:,:), intent(in)    :: sn
        real, intent(in)    :: dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2

        real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: so
        integer, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: pfl

        integer :: i,j,k
        real :: up, vp
        real :: sxplus, sxmins, syplus, symins
        real :: delta_t

        pfl = (int(sign(1.0,lmda)) + 1)/2

        so = s

        delta_t = dx/2.0d0

        !solve advection equation to get advected level set s with well defined sn  
        k = 1
        do j = jy1-NGUARD+1,jy2+NGUARD-1
          do i = ix1-NGUARD+1,ix2+NGUARD-1

          !normal vectors
          up = u(i,j,k)
          vp = v(i,j,k)

          !gradients
          sxplus = (so(i+1,j,k) - so(i,j,k)) / dx
          sxmins = (so(i,j,k) - so(i-1,j,k)) / dx
          syplus = (so(i,j+1,k) - so(i,j,k)) / dy
          symins = (so(i,j,k) - so(i,j-1,k)) / dy

          ! use dx/2 as dt to advect level set
          s(i,j,k) = so(i,j,k) + pfl(i,j,k)*delta_t*(sn(i,j,k) &
                                - max(up,0.0d0)*sxmins - min(up,0.0d0)*sxplus &
                                - max(vp,0.0d0)*symins - min(vp,0.0d0)*syplus)

          end do
        end do

        end subroutine ib_levelset_linearprojection
