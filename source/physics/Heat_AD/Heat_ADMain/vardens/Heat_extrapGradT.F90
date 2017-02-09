subroutine Heat_extrapGradT(Tnl,Tnv,T,s,pf,dx,dy,dz,nx,ny,ix1,ix2,jy1,jy2,Tnl_res,Tnv_res)

#include "Flash.h"

   implicit none
   real, dimension(:,:,:), intent(inout) :: Tnl,Tnv
   real, dimension(:,:,:), intent(in) :: T,s,pf,nx,ny
   real, intent(in) :: dx,dy,dz
   integer, intent(in) :: ix1,ix2,jy1,jy2
   real, intent(out) :: Tnl_res,Tnv_res

   integer :: i,j,k
   real, dimension(NXB,NYB) :: nxconv,nyconv,nx_plus,nx_mins,ny_plus,ny_mins
   real, dimension(NXB,NYB) :: Tlx_plus,Tlx_mins,Tly_plus,Tly_mins
   real, dimension(NXB,NYB) :: Tvx_plus,Tvx_mins,Tvy_plus,Tvy_mins

   real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: Tnl_o,Tnv_o,Tnl_i,Tnv_i

   real :: dt_ext,maxN

   integer :: l_counter, v_counter, step
   logical :: int_xp, int_xm, int_yp, int_ym


   l_counter = 0
   v_counter = 0
   k = 1

   maxN = max(maxval(abs(nx)),maxval(abs(ny)))

   !dt_ext = 0.5d0*dx*maxN
   !dt_ext = 1E-6
   dt_ext = 0.5d0/(1.0d0/dx + 1.0d0/dy)

   Tnl_i = Tnl
   Tnv_i = Tnv

   Tnl_o = Tnl

   Tnv_o = Tnv

   nxconv = sign(1.,s(ix1:ix2,jy1:jy2,k))*nx(ix1:ix2,jy1:jy2,k)
   nyconv = sign(1.,s(ix1:ix2,jy1:jy2,k))*ny(ix1:ix2,jy1:jy2,k)
   
   nx_plus = max(nxconv,0.)
   nx_mins = min(nxconv,0.)

   ny_plus = max(nyconv,0.)
   ny_mins = min(nyconv,0.)

   Tlx_plus = (Tnl_o(ix1+1:ix2+1,jy1:jy2,k)-Tnl_o(ix1:ix2,jy1:jy2,k))/dx
   Tlx_mins = (Tnl_o(ix1:ix2,jy1:jy2,k)-Tnl_o(ix1-1:ix2-1,jy1:jy2,k))/dx

   Tly_plus = (Tnl_o(ix1:ix2,jy1+1:jy2+1,k)-Tnl_o(ix1:ix2,jy1:jy2,k))/dy
   Tly_mins = (Tnl_o(ix1:ix2,jy1:jy2,k)-Tnl_o(ix1:ix2,jy1-1:jy2-1,k))/dy

   Tvx_plus = (Tnv_o(ix1+1:ix2+1,jy1:jy2,k)-Tnv_o(ix1:ix2,jy1:jy2,k))/dx
   Tvx_mins = (Tnv_o(ix1:ix2,jy1:jy2,k)-Tnv_o(ix1-1:ix2-1,jy1:jy2,k))/dx

   Tvy_plus = (Tnv_o(ix1:ix2,jy1+1:jy2+1,k)-Tnv_o(ix1:ix2,jy1:jy2,k))/dy
   Tvy_mins = (Tnv_o(ix1:ix2,jy1:jy2,k)-Tnv_o(ix1:ix2,jy1-1:jy2-1,k))/dy

   Tnl(ix1:ix2,jy1:jy2,k) = Tnl_i(ix1:ix2,jy1:jy2,k) + dt_ext*pf(ix1:ix2,jy1:jy2,k)*(-nx_mins*Tlx_plus-nx_plus*Tlx_mins &
                                                                                     -ny_mins*Tly_plus-ny_plus*Tly_mins)

   Tnv(ix1:ix2,jy1:jy2,k) = Tnv_i(ix1:ix2,jy1:jy2,k) + dt_ext*(1.0-pf(ix1:ix2,jy1:jy2,k))*(-nx_mins*Tvx_plus-nx_plus*Tvx_mins &
                                                                                           -ny_mins*Tvy_plus-ny_plus*Tvy_mins)
 
   !do j = jy1,jy2
   !  do i = ix1,ix2

   !     int_xm = .FALSE.
   !     int_xp = .FALSE.
   !     int_ym = .FALSE.
   !     int_yp = .FALSE.

   !     if(s(i,j,k)*s(i-1,j,k) .le. 0.) int_xm = .TRUE.
   !     if(s(i,j,k)*s(i+1,j,k) .le. 0.) int_xp = .TRUE.
   !     if(s(i,j,k)*s(i,j-1,k) .le. 0.) int_ym = .TRUE.
   !     if(s(i,j,k)*s(i,j+1,k) .le. 0.) int_yp = .TRUE.


   !    if(int_xm .or. int_xp .or. int_ym .or. int_yp) then

   !     Tnl(i,j,k) = Tnl_i(i,j,k) + dt_ext*pf(i,j,k)*(-nx_mins(i-ix1+1,j-jy1+1)*Tlx_plus(i-ix1+1,j-jy1+1) &
   !                                                   -nx_plus(i-ix1+1,j-jy1+1)*Tlx_mins(i-ix1+1,j-jy1+1) &
   !                                                   -ny_mins(i-ix1+1,j-jy1+1)*Tly_plus(i-ix1+1,j-jy1+1) &
   !                                                   -ny_plus(i-ix1+1,j-jy1+1)*Tly_mins(i-ix1+1,j-jy1+1))

   !     Tnv(i,j,k) = Tnv_i(i,j,k) + dt_ext*(1.0-pf(i,j,k))*(-nx_mins(i-ix1+1,j-jy1+1)*Tvx_plus(i-ix1+1,j-jy1+1) &
   !                                                         -nx_plus(i-ix1+1,j-jy1+1)*Tvx_mins(i-ix1+1,j-jy1+1) &
   !                                                         -ny_mins(i-ix1+1,j-jy1+1)*Tvy_plus(i-ix1+1,j-jy1+1) &
   !                                                         -ny_plus(i-ix1+1,j-jy1+1)*Tvy_mins(i-ix1+1,j-jy1+1))

   !    endif

   !  end do
  !end do

   Tnl_res = sum(sum(sum((Tnl_o(:,:,:)-Tnl(:,:,:))**2,1),1))
   Tnv_res = sum(sum(sum((Tnv_o(:,:,:)-Tnv(:,:,:))**2,1),1))

   Tnl_res = sqrt(Tnl_res/size(Tnl))
   Tnv_res = sqrt(Tnv_res/size(Tnv))

end subroutine Heat_extrapGradT
