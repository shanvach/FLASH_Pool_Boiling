subroutine Heat_extrapGradT(Tnl,Tnv,T,s,pf,dx,dy,dz,nx,ny,ix1,ix2,jy1,jy2,Tnl_res,Tnv_res)

#include "Flash.h"

    implicit none
    real, dimension(:,:,:), intent(inout) :: Tnl,Tnv
    real, dimension(:,:,:), intent(in) :: T,s,pf,nx,ny
    real, intent(in) :: dx,dy,dz
    integer, intent(in) :: ix1,ix2,jy1,jy2
    real, intent(out) :: Tnl_res,Tnv_res

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

    !dt_ext = 0.5d0*dx*maxN
    !dt_ext = 1E-6
    dt_ext = 0.5d0/(1.0d0/dx + 1.0d0/dy)

    Tnl_i = Tnl
    Tnv_i = Tnv

    Tnl_o = Tnl

    Tnv_o = Tnv

    k = 1
    do i=ix1,ix2
     do j=jy1,jy2

      if(pf(i,j,k) .eq. 1.) then

          nxconv = nx(i,j,k)
          nyconv = ny(i,j,k)

          nx_plus = max(nxconv,0.)
          nx_mins = min(nxconv,0.)

          ny_plus = max(nyconv,0.)
          ny_mins = min(nyconv,0.)

          Tlx_plus = (Tnl_o(i+1,j,k) - Tnl_o(i,j,k))/dx
          Tlx_mins = (Tnl_o(i,j,k) - Tnl_o(i-1,j,k))/dx

          Tly_plus = (Tnl_o(i,j+1,k) - Tnl_o(i,j,k))/dy
          Tly_mins = (Tnl_o(i,j,k) - Tnl_o(i,j-1,k))/dy

          Tnl(i,j,k) = Tnl_i(i,j,k) + dt_ext*(-nx_mins*Tlx_plus-nx_plus*Tlx_mins &
                                                    -ny_mins*Tly_plus-ny_plus*Tly_mins)

      else if(pf(i,j,k) .eq. 0.) then

          nxconv = -nx(i,j,k)
          nyconv = -ny(i,j,k)

          nx_plus = max(nxconv,0.)
          nx_mins = min(nxconv,0.)

          ny_plus = max(nyconv,0.)
          ny_mins = min(nyconv,0.)

          Tvx_plus = (Tnv_o(i+1,j,k) - Tnv_o(i,j,k))/dx
          Tvx_mins = (Tnv_o(i,j,k) - Tnv_o(i-1,j,k))/dx

          Tvy_plus = (Tnv_o(i,j+1,k) - Tnv_o(i,j,k))/dy
          Tvy_mins = (Tnv_o(i,j,k) - Tnv_o(i,j-1,k))/dy

          Tnv(i,j,k) = Tnv_i(i,j,k) + dt_ext*(-nx_mins*Tvx_plus-nx_plus*Tvx_mins &
                                                    -ny_mins*Tvy_plus-ny_plus*Tvy_mins)

       end if

          !nx_plus = max(nxconv,0.)
          !nx_mins = min(nxconv,0.)

          !ny_plus = max(nyconv,0.)
          !ny_mins = min(nyconv,0.)

          !Tlx_plus = (Tnl_o(i+1,j,k) - Tnl_o(i,j,k))/dx
          !Tlx_mins = (Tnl_o(i,j,k) - Tnl_o(i-1,j,k))/dx

          !Tly_plus = (Tnl_o(i,j+1,k) - Tnl_o(i,j,k))/dy
          !Tly_mins = (Tnl_o(i,j,k) - Tnl_o(i,j-1,k))/dy

          !Tvx_plus = (Tnv_o(i+1,j,k) - Tnv_o(i,j,k))/dx
          !Tvx_mins = (Tnv_o(i,j,k) - Tnv_o(i-1,j,k))/dx

          !Tvy_plus = (Tnv_o(i,j+1,k) - Tnv_o(i,j,k))/dy
          !Tvy_mins = (Tnv_o(i,j,k) - Tnv_o(i,j-1,k))/dy

          !Tnl(i,j,k) = Tnl_i(i,j,k) + dt_ext*(-nx_mins*Tlx_plus-nx_plus*Tlx_mins&
          !                                          -ny_mins*Tly_plus-ny_plus*Tly_mins)

          !Tnv(i,j,k) = Tnv_i(i,j,k) + dt_ext*(-nx_mins*Tvx_plus-nx_plus*Tvx_mins&
          !                                          -ny_mins*Tvy_plus-ny_plus*Tvy_mins)

     end do
   end do
  
   Tnl_res = sum(sum(sum((Tnl_o(:,:,:)-Tnl(:,:,:))**2,1),1))
   Tnv_res = sum(sum(sum((Tnv_o(:,:,:)-Tnv(:,:,:))**2,1),1))

   Tnl_res = sqrt(Tnl_res/size(Tnl))
   Tnv_res = sqrt(Tnv_res/size(Tnv))

end subroutine Heat_extrapGradT
