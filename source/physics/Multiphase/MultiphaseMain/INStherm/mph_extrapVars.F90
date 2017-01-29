subroutine mph_extrapVars(u,v,T,p,s,pf,dx,dy,dz,nx,ny,ix1,ix2,jy1,jy2)

#include "Flash.h"

    implicit none
    real, dimension(:,:,:), intent(inout) :: u,v
    real, dimension(:,:,:), intent(in) :: T,p,s,pf,nx,ny
    real, intent(in) :: dx,dy,dz
    integer, intent(in) :: ix1,ix2,jy1,jy2

    integer :: i,j,k
    real :: nxconv,nyconv,nx_plus,nx_mins,ny_plus,ny_mins
    real :: Tlx_plus, Tlx_mins, Tly_plus, Tly_mins
    real :: Tvx_plus, Tvx_mins, Tvy_plus, Tvy_mins

    real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: Tnl_o,Tnv_o,Tnl_i,Tnv_i

    real :: dt_ext,maxN

    integer :: l_counter, v_counter, step

    l_counter = 0
    v_counter = 0

    maxN = max(maxval(abs(nx)),maxval(abs(ny)))

    dt_ext = 0.5d0*dx*maxN
    !dt_ext = 1E-6
    !dt_ext = 0.5d0/(1.0d0/dx + 1.0d0/dy)

end subroutine mph_extrapVars
