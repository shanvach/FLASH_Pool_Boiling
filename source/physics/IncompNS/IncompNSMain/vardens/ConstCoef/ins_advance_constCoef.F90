subroutine ins_rhsConstCoef (ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz,dt,rho1,rho2,p,&
                             po,sigp,sigox,sigoy,sigoz,sigoox,sigooy,sigooz,rho1x,rho2x,&
                             rho1y,rho2y,rho1z,rho2z,rhoz,pf,s,divu,blockID)

#include "Flash.h"
#include "constants.h"

        implicit none
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2,blockID
        real, intent(in) :: dt,rho1,rho2,rhoz

        real, dimension(:,:), intent(in) :: dx,dy,dz

        real, dimension(:,:,:), intent(in) :: p,po,sigp,sigox,sigoy,sigoz,&
                                              sigoox,sigooy,sigooz,pf,s

        real, dimension(:,:,:), intent(inout) :: divu,rho1x,rho2x,&
                                                 rho1y,rho2y,rho1z,rho2z
                                                 
        integer :: i,j,k

        real :: d1x,d1y,p1x1,p1x2,p1x,p1y1,p1y2,p1y,&
                d2x,d2y,p2x1,p2x2,p2x,p2y1,p2y2,p2y,&
                d1z,d2z,p1z1,p1z2,p1z,p2z1,p2z2,p2z,&
                term1,term2,term3
        
        real :: th,adensx,adensy,adensz

        do k=kz1,kz2
        do i=ix1,ix2
           do j=jy1,jy2

              term1 = sigp(i,j,k)
        
              adensx = rho1x(i+1,j,k)+rho2x(i+1,j,k)
              
              d1x = 1.0-(rhoz*adensx)
              p1x =((2.0*p(i+1,j,k)-1.0*po(i+1,j,k))-(2.0*p(i,j,k)-1.0*po(i,j,k)))*dx(i,RIGHT_EDGE) -&
                    (2*sigox(i+1,j,k)-sigoox(i+1,j,k))

              adensx = rho1x(i,j,k)+rho2x(i,j,k)

              d2x = 1.0-(rhoz*adensx)
              p2x =((2.0*p(i,j,k)-1.0*po(i,j,k))-(2.0*p(i-1,j,k)-1.0*po(i-1,j,k)))*dx(i,LEFT_EDGE)-&
                    (2*sigox(i,j,k)-sigoox(i,j,k))

              adensy = rho1y(i,j+1,k)+rho2y(i,j+1,k)

              d1y = 1.0-(rhoz*adensy)
              p1y =((2.0*p(i,j+1,k)-1.0*po(i,j+1,k))-(2.0*p(i,j,k)-1.0*po(i,j,k)))*dy(j,RIGHT_EDGE)-&
                    (2*sigoy(i,j+1,k)-sigooy(i,j+1,k))

              adensy = rho1y(i,j,k)+rho2y(i,j,k)

              d2y = 1.0-(rhoz*adensy)
              p2y =((2.0*p(i,j,k)-1.0*po(i,j,k))-(2.0*p(i,j-1,k)-1.0*po(i,j-1,k)))*dy(j,LEFT_EDGE)-&
                    (2*sigoy(i,j,k)-sigooy(i,j,k))

              term2 = (((d1x*p1x)-(d2x*p2x))*dx(i,CENTER))+(((d1y*p1y)-(d2y*p2y))*dy(j,CENTER))

#if NDIM == 3
              adensz = rho1z(i,j,k+1)+rho2z(i,j,k+1)

              d1z = 1.0-(rhoz*adensz)
              p1z =((2.0*p(i,j,k+1)-1.0*po(i,j,k+1))-(2.0*p(i,j,k)-1.0*po(i,j,k)))*dz(k,RIGHT_EDGE)-&
                    (2*sigoz(i,j,k+1)-sigooz(i,j,k+1))

              adensz = rho1z(i,j,k)+rho2z(i,j,k)

              d2z = 1.0-(rhoz*adensz)
              p2z =((2.0*p(i,j,k)-1.0*po(i,j,k))-(2.0*p(i,j,k-1)-1.0*po(i,j,k-1)))*dz(k,LEFT_EDGE)-&
                    (2*sigoz(i,j,k)-sigooz(i,j,k))
 
              term2 = term2 + (((d1z*p1z)-(d2z*p2z))*dz(k,CENTER))
#endif

              term3 = (rhoz/dt)*divu(i,j,k)

              divu(i,j,k) = term1+term2+term3
 
          enddo
        enddo
        enddo

end subroutine ins_rhsConstCoef

subroutine ins_correctorConstCoef(uni,vni,wni,sigx,sigy,sigz,sigox,sigoy,sigoz,&
                                  sigoox,sigooy,sigooz,rho1x,rho2x,rho1y,rho2y,&
                                  rho1z,rho2z,p,po,poo,ix1,ix2,jy1,jy2,kz1,kz2,&
                                  rhoz,dt,dx,dy,dz,rho1,rho2,pf,s,alfa,blockID)
#include "Flash.h"

        implicit none

        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2,blockID
        real, intent(in) :: dt,alfa,rhoz,rho1,rho2

        real, dimension(:), intent(in) :: dx,dy,dz

        real, dimension(:,:,:), intent(in) :: p,po,poo,sigx,sigy,sigz,sigox,sigoy,sigoz,sigoox,&
                                              sigooy,sigooz,pf,s

        real, dimension(:,:,:), intent(inout) :: uni,vni,wni,rho1x,rho2x,&
                                                 rho1y,rho2y,rho1z,rho2z

        real :: Mdens,coef,coef1,coef2,coef3,irhoz,th
        integer :: i,j,k

        irhoz = 1.0/rhoz
       
        do k=kz1,kz2
           do j=jy1,jy2
              do i=ix1,ix2+1

              Mdens = (rho1x(i,j,k) + rho2x(i,j,k))

              uni(i,j,k) = uni(i,j,k) - dt*&
                           (irhoz*(((p(i,j,k)-p(i-1,j,k))*dx(i)) - sigx(i,j,k) )+&
                           (Mdens-irhoz)*(( (2.0*po(i,j,k)-2.0*po(i-1,j,k))*dx(i)-&
                                             (1.0*poo(i,j,k)-1.0*poo(i-1,j,k))*dx(i) ) -&
                                             (2.0*sigox(i,j,k)-1.0*sigoox(i,j,k)) ))

                enddo
            enddo
        enddo
 
         do k=kz1,kz2
           do j=jy1,jy2+1
              do i=ix1,ix2

              Mdens = (rho1y(i,j,k) + rho2y(i,j,k))

              vni(i,j,k) = vni(i,j,k) -  dt*&
                           (irhoz*(((p(i,j,k)-p(i,j-1,k))*dy(j)) - sigy(i,j,k) )+&
                           (Mdens-irhoz)*(( (2.0*po(i,j,k)-2.0*po(i,j-1,k))*dy(j)-&
                                             (1.0*poo(i,j,k)-1.0*poo(i,j-1,k))*dy(j) ) -&
                                             (2.0*sigoy(i,j,k)-1.0*sigooy(i,j,k)) ))

                enddo
            enddo
        enddo

#if NDIM == 3
        do k=kz1,kz2+1
           do j=jy1,jy2
              do i=ix1,ix2

              Mdens = (rho1z(i,j,k) + rho2z(i,j,k))

              wni(i,j,k) = wni(i,j,k) - dt*&
                           (irhoz*(((p(i,j,k)-p(i,j,k-1))*dz(k)) - sigz(i,j,k) )+&
                           (Mdens-irhoz)*(( (2.0*po(i,j,k)-2.0*po(i,j,k-1))*dz(k)-&
                                             (1.0*poo(i,j,k)-1.0*poo(i,j,k-1))*dz(k) ) -&
                                             (2.0*sigoz(i,j,k)-1.0*sigooz(i,j,k)) ))

                enddo
            enddo
        enddo
#endif

end subroutine ins_correctorConstCoef

subroutine ins_correctorConstCoef_IB(uni,vni,wni,sigx,sigy,sigz,sigox,sigoy,sigoz,&
                                  sigoox,sigooy,sigooz,rho1x,rho2x,rho1y,rho2y,&
                                  rho1z,rho2z,p,po,poo,ix1,ix2,jy1,jy2,kz1,kz2,&
                                  rhoz,dt,dx,dy,dz,rho1,rho2,pf,s,alfa,&
                                  poldx,poldy,poldz,blockID)
#include "Flash.h"

        implicit none

        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2,blockID
        real, intent(in) :: dt,alfa,rhoz,rho1,rho2

        real, dimension(:), intent(in) :: dx,dy,dz

        real, dimension(:,:,:), intent(in) :: p,po,poo,sigx,sigy,sigz,sigox,sigoy,sigoz,sigoox,&
                                              sigooy,sigooz,pf,s

        real, dimension(:,:,:), intent(inout) :: uni,vni,wni,rho1x,rho2x,&
                                                 rho1y,rho2y,rho1z,rho2z,&
                                                 poldx,poldy,poldz

        real :: Mdens,coef,coef1,coef2,coef3,irhoz,th
        integer :: i,j,k

        irhoz = 1.0/rhoz

        poldx = 0.0
        poldy = 0.0
        poldz = 0.0
 
        do k=kz1,kz2
           do j=jy1,jy2
              do i=ix1,ix2+1

              Mdens = (rho1x(i,j,k) + rho2x(i,j,k))

              uni(i,j,k) = uni(i,j,k) - dt*alfa*&
                           (irhoz*(((p(i,j,k)-p(i-1,j,k))*dx(i)) - sigx(i,j,k) )+&
                           (Mdens-irhoz)*(( (2.0*po(i,j,k)-2.0*po(i-1,j,k))*dx(i)-&
                                             (1.0*poo(i,j,k)-1.0*poo(i-1,j,k))*dx(i) ) -&
                                             (2.0*sigox(i,j,k)-1.0*sigoox(i,j,k)) ))
                 
              poldx(i,j,k) = irhoz*(((p(i,j,k)-p(i-1,j,k))*dx(i)) - sigx(i,j,k) )+&
                            (Mdens-irhoz)*(( (2.0*po(i,j,k)-2.0*po(i-1,j,k))*dx(i)-&
                                             (1.0*poo(i,j,k)-1.0*poo(i-1,j,k))*dx(i) ) -&
                                             (2.0*sigox(i,j,k)-1.0*sigoox(i,j,k)) )

                enddo
            enddo
        enddo
 
         do k=kz1,kz2
           do j=jy1,jy2+1
              do i=ix1,ix2

              Mdens = (rho1y(i,j,k) + rho2y(i,j,k))

              vni(i,j,k) = vni(i,j,k) -  dt*alfa*&
                           (irhoz*(((p(i,j,k)-p(i,j-1,k))*dy(j)) - sigy(i,j,k) )+&
                           (Mdens-irhoz)*(( (2.0*po(i,j,k)-2.0*po(i,j-1,k))*dy(j)-&
                                             (1.0*poo(i,j,k)-1.0*poo(i,j-1,k))*dy(j) ) -&
                                             (2.0*sigoy(i,j,k)-1.0*sigooy(i,j,k)) ))
    
              poldy(i,j,k) = irhoz*(((p(i,j,k)-p(i,j-1,k))*dy(j)) - sigy(i,j,k) )+&
                            (Mdens-irhoz)*(( (2.0*po(i,j,k)-2.0*po(i,j-1,k))*dy(j)-&
                                             (1.0*poo(i,j,k)-1.0*poo(i,j-1,k))*dy(j) ) -&
                                             (2.0*sigoy(i,j,k)-1.0*sigooy(i,j,k)) )

                enddo
            enddo
        enddo

#if NDIM == 3
        do k=kz1,kz2+1
           do j=jy1,jy2
              do i=ix1,ix2

              Mdens = (rho1z(i,j,k) + rho2z(i,j,k))

              wni(i,j,k) = wni(i,j,k) - dt*alfa*&
                           (irhoz*(((p(i,j,k)-p(i,j,k-1))*dz(k)) - sigz(i,j,k) )+&
                           (Mdens-irhoz)*(( (2.0*po(i,j,k)-2.0*po(i,j,k-1))*dz(k)-&
                                             (1.0*poo(i,j,k)-1.0*poo(i,j,k-1))*dz(k) ) -&
                                             (2.0*sigoz(i,j,k)-1.0*sigooz(i,j,k)) ))
              
              poldz(i,j,k) = irhoz*(((p(i,j,k)-p(i,j,k-1))*dz(k)) - sigz(i,j,k) )+&
                            (Mdens-irhoz)*(( (2.0*po(i,j,k)-2.0*po(i,j,k-1))*dz(k)-&
                                             (1.0*poo(i,j,k)-1.0*poo(i,j,k-1))*dz(k) ) -&
                                             (2.0*sigoz(i,j,k)-1.0*sigooz(i,j,k)) )
 
                enddo
            enddo
        enddo
#endif

end subroutine ins_correctorConstCoef_IB
