subroutine Heat_getQin(u,v,pf,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2,qin)


#include "Flash.h"

        implicit none

        real,intent(in),dimension(:,:,:) :: u,v,pf
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
                           if(pf(i,j,k) .eq. 1 .and. pf(i+1,j,k) .eq. 0) qin = qin + u(i+1,j,k)*dxdz
                           if(pf(i,j,k) .eq. 0 .and. pf(i+1,j,k) .eq. 1) qin = qin - u(i+1,j,k)*dxdz
                           if(pf(i,j,k) .eq. 1 .and. pf(i,j+1,k) .eq. 0) qin = qin + v(i,j+1,k)*dydz
                           if(pf(i,j,k) .eq. 0 .and. pf(i,j+1,k) .eq. 1) qin = qin - v(i,j+1,k)*dydz
                        end do
                end do
        end do

end subroutine Heat_getQin
