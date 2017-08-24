subroutine Heat_extrapGradT(Tnl,Tnv,T,s,pf,dx,dy,dz,nx,ny,ix1,ix2,jy1,jy2,Tnl_res,Tnv_res,mflg)

#include "Flash.h"

#define FIRST_ORDER_UPWIND
!#define SECOND_ORDER_UPWIND
!#define THIRD_ORDER_UPWIND

   implicit none
   real, dimension(:,:,:), intent(inout) :: Tnl,Tnv
   real, dimension(:,:,:), intent(in) :: T,s,pf,nx,ny,mflg
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

   dt_ext = 0.5d0*(dx/maxN)
   !dt_ext = 1E-6
   !dt_ext = 0.5d0/(1.0d0/dx + 1.0d0/dy)
   !dt_ext = dx

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

#ifdef FIRST_ORDER_UPWIND
   !---First Order Upwind-------!
   Tlx_plus = (Tnl_o(ix1+1:ix2+1,jy1:jy2,k)-Tnl_o(ix1:ix2,jy1:jy2,k))/dx
   Tlx_mins = (Tnl_o(ix1:ix2,jy1:jy2,k)-Tnl_o(ix1-1:ix2-1,jy1:jy2,k))/dx

   Tly_plus = (Tnl_o(ix1:ix2,jy1+1:jy2+1,k)-Tnl_o(ix1:ix2,jy1:jy2,k))/dy
   Tly_mins = (Tnl_o(ix1:ix2,jy1:jy2,k)-Tnl_o(ix1:ix2,jy1-1:jy2-1,k))/dy

   Tvx_plus = (Tnv_o(ix1+1:ix2+1,jy1:jy2,k)-Tnv_o(ix1:ix2,jy1:jy2,k))/dx
   Tvx_mins = (Tnv_o(ix1:ix2,jy1:jy2,k)-Tnv_o(ix1-1:ix2-1,jy1:jy2,k))/dx

   Tvy_plus = (Tnv_o(ix1:ix2,jy1+1:jy2+1,k)-Tnv_o(ix1:ix2,jy1:jy2,k))/dy
   Tvy_mins = (Tnv_o(ix1:ix2,jy1:jy2,k)-Tnv_o(ix1:ix2,jy1-1:jy2-1,k))/dy
   !-----------------------------!
#endif

#ifdef SECOND_ORDER_UPWIND
   !---Second Order Upwind-------!
   Tlx_plus = (+4*Tnl_o(ix1+1:ix2+1,jy1:jy2,k)-3*Tnl_o(ix1:ix2,jy1:jy2,k)-Tnl_o(ix1+2:ix2+2,jy1:jy2,k))/(2.0*dx)
   Tlx_mins = (-4*Tnl_o(ix1-1:ix2-1,jy1:jy2,k)+3*Tnl_o(ix1:ix2,jy1:jy2,k)+Tnl_o(ix1-2:ix2-2,jy1:jy2,k))/(2.0*dx)

   Tly_plus = (+4*Tnl_o(ix1:ix2,jy1+1:jy2+1,k)-3*Tnl_o(ix1:ix2,jy1:jy2,k)-Tnl_o(ix1:ix2,jy1+2:jy2+2,k))/(2.0*dy)
   Tly_mins = (-4*Tnl_o(ix1:ix2,jy1-1:jy2-1,k)+3*Tnl_o(ix1:ix2,jy1:jy2,k)+Tnl_o(ix1:ix2,jy1-2:jy2-2,k))/(2.0*dy)

   Tvx_plus = (+4*Tnv_o(ix1+1:ix2+1,jy1:jy2,k)-3*Tnv_o(ix1:ix2,jy1:jy2,k)-Tnv_o(ix1+2:ix2+2,jy1:jy2,k))/(2.0*dx)
   Tvx_mins = (-4*Tnv_o(ix1-1:ix2-1,jy1:jy2,k)+3*Tnv_o(ix1:ix2,jy1:jy2,k)+Tnv_o(ix1-2:ix2-2,jy1:jy2,k))/(2.0*dx)

   Tvy_plus = (+4*Tnv_o(ix1:ix2,jy1+1:jy2+1,k)-3*Tnv_o(ix1:ix2,jy1:jy2,k)-Tnv_o(ix1:ix2,jy1+2:jy2+2,k))/(2.0*dy)
   Tvy_mins = (-4*Tnv_o(ix1:ix2,jy1-1:jy2-1,k)+3*Tnv_o(ix1:ix2,jy1:jy2,k)+Tnv_o(ix1:ix2,jy1-2:jy2-2,k))/(2.0*dy)
   !-----------------------------!
#endif

#ifdef THIRD_ORDER_UPWIND
   !---Third Order Upwind--------!
   Tlx_plus = (-2*Tnl_o(ix1-1:ix2-1,jy1:jy2,k)+6*Tnl_o(ix1+1:ix2+1,jy1:jy2,k) &
               -3*Tnl_o(ix1:ix2,jy1:jy2,k)-Tnl_o(ix1+2:ix2+2,jy1:jy2,k))/(6.0*dx)

   Tlx_mins = (+2*Tnl_o(ix1+1:ix2+1,jy1:jy2,k)-6*Tnl_o(ix1-1:ix2-1,jy1:jy2,k) &
               +3*Tnl_o(ix1:ix2,jy1:jy2,k)+Tnl_o(ix1-2:ix2-2,jy1:jy2,k))/(6.0*dx)

   Tly_plus = (-2*Tnl_o(ix1:ix2,jy1-1:jy2-1,k)+6*Tnl_o(ix1:ix2,jy1+1:jy2+1,k) &
               -3*Tnl_o(ix1:ix2,jy1:jy2,k)-Tnl_o(ix1:ix2,jy1+2:jy2+2,k))/(6.0*dy)

   Tly_mins = (+2*Tnl_o(ix1:ix2,jy1+1:jy2+1,k)-6*Tnl_o(ix1:ix2,jy1-1:jy2-1,k) &
               +3*Tnl_o(ix1:ix2,jy1:jy2,k)+Tnl_o(ix1:ix2,jy1-2:jy2-2,k))/(6.0*dy)

   Tvx_plus = (-2*Tnv_o(ix1-1:ix2-1,jy1:jy2,k)+6*Tnv_o(ix1+1:ix2+1,jy1:jy2,k) &
               -3*Tnv_o(ix1:ix2,jy1:jy2,k)-Tnv_o(ix1+2:ix2+2,jy1:jy2,k))/(6.0*dx)

   Tvx_mins = (+2*Tnv_o(ix1+1:ix2+1,jy1:jy2,k)-6*Tnv_o(ix1-1:ix2-1,jy1:jy2,k) &
               +3*Tnv_o(ix1:ix2,jy1:jy2,k)+Tnv_o(ix1-2:ix2-2,jy1:jy2,k))/(6.0*dx)

   Tvy_plus = (-2*Tnv_o(ix1:ix2,jy1-1:jy2-1,k)+6*Tnv_o(ix1:ix2,jy1+1:jy2+1,k) &
               -3*Tnv_o(ix1:ix2,jy1:jy2,k)-Tnv_o(ix1:ix2,jy1+2:jy2+2,k))/(6.0*dy)

   Tvy_mins = (+2*Tnv_o(ix1:ix2,jy1+1:jy2+1,k)-6*Tnv_o(ix1:ix2,jy1-1:jy2-1,k) &
               +3*Tnv_o(ix1:ix2,jy1:jy2,k)+Tnv_o(ix1:ix2,jy1-2:jy2-2,k))/(6.0*dy)
   !-----------------------------!
#endif

   Tnl(ix1:ix2,jy1:jy2,k) = Tnl_i(ix1:ix2,jy1:jy2,k) + dt_ext*pf(ix1:ix2,jy1:jy2,k)*(-nx_mins*Tlx_plus-nx_plus*Tlx_mins &
                                                                                     -ny_mins*Tly_plus-ny_plus*Tly_mins)

   Tnv(ix1:ix2,jy1:jy2,k) = Tnv_i(ix1:ix2,jy1:jy2,k) + dt_ext*(1.0-pf(ix1:ix2,jy1:jy2,k))*(-nx_mins*Tvx_plus-nx_plus*Tvx_mins &
                                                                                           -ny_mins*Tvy_plus-ny_plus*Tvy_mins)

   Tnl_res = sum(sum(sum((mflg*(Tnl_o(:,:,:)-Tnl(:,:,:)))**2,1),1))
   Tnv_res = sum(sum(sum((mflg*(Tnv_o(:,:,:)-Tnv(:,:,:)))**2,1),1))

   Tnl_res = sqrt(Tnl_res/size(Tnl))
   Tnv_res = sqrt(Tnv_res/size(Tnv))

end subroutine Heat_extrapGradT
