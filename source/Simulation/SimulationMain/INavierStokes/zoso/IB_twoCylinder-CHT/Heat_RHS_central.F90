subroutine Heat_RHS_central(T_rhs, T_o, T_g, u, v, dx, dy, dz,inRe, ix1,ix2, jy1,jy2,&
                    rho1x,rho2x,rho1y,rho2y,alph,pf,s,mdot,nrmx,nrmy,smrh,curv,lambda,phip,tfrx,tfry)

  use Heat_AD_data
  use Multiphase_data, only: mph_cp2,mph_thco2, mph_rho2,mph_rho1
  use Heat_AD_interface, only: Heat_GFMstencil_o1, Heat_GFMstencil_o2
  use Driver_data, only: dr_dt
 
#include "Heat_AD.h"

  implicit none
  real, dimension(:,:,:), intent(inout) :: T_rhs
  real, dimension(:,:,:), intent(in) :: T_o, T_g, tfrx, tfry
  real, dimension(:,:,:), intent(in) :: u,v
  real, intent(in) :: dx, dy, dz, inRe
  integer, intent(in) :: ix1, ix2, jy1, jy2
  real, dimension(:,:,:),intent(in) :: rho1x,rho2x,rho1y,rho2y,alph
  real, dimension(:,:,:),intent(in) :: pf,s,mdot,nrmx,nrmy,smrh,curv,lambda,phip

  real :: T_res

  integer :: i,j,k

  real :: Tx_plus, Tx_mins, Ty_plus, Ty_mins, Tij
  real :: ul, ur, vl, vr, uc, vc
  real :: Txx, Tyy
  real :: coeff
  real :: tol
  real :: rhoc,rhoxm,rhoym,rhoxp,rhoyp
  real :: thxp1, thxm1
  real :: thyp1, thym1
  real :: thxp2, thxm2
  real :: thyp2, thym2

  real :: eps, &
          s1r,s2r,s3r,s4r,s5r,s1l,s2l,s3l,s4l,s5l, &
          rIS1r,rIS2r,rIS3r,rIS1l,rIS2l,rIS3l, &
          aT1r,aT2r,aT3r,aT1l,aT2l,aT3l, &
          a1r,a2r,a3r,a1l,a2l,a3l, &
          fT1r,fT2r,fT3r,fT1l,fT2l,fT3l, &
          frx,flx,fry,fly

  tol = 0.01
  eps = 1E-15

  k = 1

  do j=jy1,jy2
     do i=ix1,ix2

     ul = u(i,j,k)
     ur = u(i+1,j,k)
     vl = v(i,j,k)
     vr = v(i,j+1,k)

     uc = (ul + ur)/2.0
     vc = (vl + vr)/2.0

     coeff = inRe/ht_Pr

     Tx_plus = T_o(i+1,j,k)
     Tx_mins = T_o(i-1,j,k)

     Ty_plus = T_o(i,j+1,k)
     Ty_mins = T_o(i,j-1,k)

     Tij = T_o(i,j,k)

     thxp1 = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i+1,j,k)))
     thxm1 = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i-1,j,k)))
     thyp1 = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j+1,k)))
     thym1 = abs(s(i,j,k))/(abs(s(i,j,k))+abs(s(i,j-1,k)))

     !______________________Diffusion Terms_______________________!

     ! Case 1 !
     if(s(i,j,k)*s(i+1,j,k) .le. 0.d0 .and. lambda(i,j,k) .lt. 0.0 .and. sign(1.,s(i,j,k)) .ne. sign(1.,phip(i+1,j,k)))&
     call Heat_GFMstencil_o1(Tx_plus,Tij,ht_Tsat,max(tol,thxp1))
     ! End of Case 1 !

     ! Case 2 !
     if(s(i,j,k)*s(i-1,j,k) .le. 0.d0 .and. lambda(i,j,k) .lt. 0.0 .and. sign(1.,s(i,j,k)) .ne. sign(1.,phip(i-1,j,k))) &
     call Heat_GFMstencil_o1(Tx_mins,Tij,ht_Tsat,max(tol,thxm1))
     ! End of Case 2 !

     ! Case 3 !
     if(s(i,j,k)*s(i,j+1,k) .le. 0.d0 .and. lambda(i,j,k) .lt. 0.0 .and. sign(1.,s(i,j,k)) .ne. sign(1.,phip(i,j+1,k))) &
     call Heat_GFMstencil_o1(Ty_plus,Tij,ht_Tsat,max(tol,thyp1))
     ! End of Case 3 !

     ! Case 4 !
     if(s(i,j,k)*s(i,j-1,k) .le. 0.d0 .and. lambda(i,j,k) .lt. 0.0 .and. sign(1.,s(i,j,k)) .ne. sign(1.,phip(i,j-1,k))) &
     call Heat_GFMstencil_o1(Ty_mins,Tij,ht_Tsat,max(tol,thym1))
     ! End of Case 4 ! 
 
     !______________________IB Terms_______________________!

     !! Case 1 !
     if((lambda(i,j,k)*lambda(i+1,j,k) .le. 0.d0 .and. lambda(i,j,k) .lt. 0.0 .and. sign(1.,s(i,j,k)) .eq. sign(1.,phip(i+1,j,k))) .or. &
        (lambda(i,j,k)*lambda(i+1,j,k) .le. 0.d0 .and. lambda(i,j,k) .ge. 0.0 .and. sign(1.,s(i+1,j,k)) .eq. sign(1.,phip(i,j,k)))  )   &
     Tx_plus = T_g(i+1,j,k)
     !! End of Case 1 !

     !! Case 2 !
     if((lambda(i,j,k)*lambda(i-1,j,k) .le. 0.d0 .and. lambda(i,j,k) .lt. 0.0 .and. sign(1.,s(i,j,k)) .eq. sign(1.,phip(i-1,j,k))) .or. &
        (lambda(i,j,k)*lambda(i-1,j,k) .le. 0.d0 .and. lambda(i,j,k) .ge. 0.0 .and. sign(1.,s(i-1,j,k)) .eq. sign(1.,phip(i,j,k)))  )   &
     Tx_mins = T_g(i-1,j,k)
     !! End of Case 2 !

     !! Case 3 !
     if((lambda(i,j,k)*lambda(i,j+1,k) .le. 0.d0 .and. lambda(i,j,k) .lt. 0.0 .and. sign(1.,s(i,j,k)) .eq. sign(1.,phip(i,j+1,k))) .or. &
        (lambda(i,j,k)*lambda(i,j+1,k) .le. 0.d0 .and. lambda(i,j,k) .ge. 0.0 .and. sign(1.,s(i,j+1,k)) .eq. sign(1.,phip(i,j,k)))  )   &
     Ty_plus = T_g(i,j+1,k)
     !! End of Case 3 !

     !! Case 4 !
     if((lambda(i,j,k)*lambda(i,j-1,k) .le. 0.d0 .and. lambda(i,j,k) .lt. 0.0 .and. sign(1.,s(i,j,k)) .eq. sign(1.,phip(i,j-1,k))) .or. &
        (lambda(i,j,k)*lambda(i,j-1,k) .le. 0.d0 .and. lambda(i,j,k) .ge. 0.0 .and. sign(1.,s(i,j-1,k)) .eq. sign(1.,phip(i,j,k)))  )   &
     Ty_mins = T_g(i,j-1,k)
     !! End of Case 4 ! 
    
    !_______________________________RHS TERM______________________________________!

    Txx = alph(i,j,k)*(coeff*(Tx_plus-Tij)/dx - coeff*(Tij-Tx_mins)/dx)/dx
    Tyy = alph(i,j,k)*(coeff*(Ty_plus-Tij)/dy - coeff*(Tij-Ty_mins)/dy)/dy

    if(s(i,j,k) .ge. 0.0) then

     T_rhs(i,j,k) = 0.0

    else

        if(lambda(i,j,k) .lt. 0.0) then

        T_rhs(i,j,k) = - uc*(Tx_plus - Tx_mins)/(2.0*dx) &
                       - vc*(Ty_plus - Ty_mins)/(2.0*dy) &
                       - uc*tfrx(i,j,k) - vc*tfry(i,j,k) &
                       + Txx + Tyy

        else

        T_rhs(i,j,k) = Txx + Tyy

        end if

    end if

    end do
  end do 

end subroutine Heat_RHS_central
