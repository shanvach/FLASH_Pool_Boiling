subroutine Heat_getQin(mdot,smrh,dfun,nrmx,nrmy,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2,qin)


#include "Flash.h"

        implicit none

        real,intent(in),dimension(:,:,:) :: mdot,smrh,dfun,nrmx,nrmy
        real,intent(in) :: dx,dy,dz
        real,intent(inout) :: qin
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2

        integer :: i,j,k

        real :: mdoti, mdotj, mdotk, densx, densy, densz, dxdy, dxdz, dydz

        dxdy = dx*dy
        dxdz = dx*dz
        dydz = dy*dz

        do k=kz1,kz2
                do j=jy1,jy2
                        do i=ix1,ix2

                                mdoti = 0.5*(mdot(i,j,k)+mdot(i+1,j,k))
                                mdotj = 0.5*(mdot(i,j,k)+mdot(i,j+1,k))
        
                                densx = 0.5*(smrh(i,j,k)+smrh(i+1,j,k))
                                densy = 0.5*(smrh(i,j,k)+smrh(i,j+1,k))                                

#if NDIM == 3
                                mdotk = 0.5*(mdot(i,j,k)+mdot(i,j,k+1))
                                densz = 0.5*(smrh(i,j,k)+smrh(i,j,k+1))
#endif        

                                qin = qin + mdot(i,j,k)*smrh(i,j,k)*(dx*nrmx(i,j,k) + dy*nrmy(i,j,k))

                                !if(dfun(i,j,k)*dfun(i+1,j,k) .le. 0.0) qin = qin - mdoti*densx*dxdz
                                !if(dfun(i,j,k)*dfun(i,j+1,k) .le. 0.0) qin = qin - mdotj*densy*dydz

#if NDIM == 3
                                !if(dfun(i,j,k)*dfun(i,j,k+1) .le. 0.0) qin = qin - mdotk*densz*dxdy
#endif


                        end do
                end do
        end do

end subroutine Heat_getQin
