subroutine Heat_getConvel(mdot,smrh,xnorm,ynorm,znorm,uconv,ix1,ix2,jy1,jy2,kz1,kz2)


#include "Flash.h"

        implicit none

        real,intent(in),dimension(:,:,:) :: mdot,smrh
        real,intent(in),dimension(:,:,:) :: xnorm,ynorm,znorm
        real,intent(inout) :: uconv
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2

        integer :: i,j,k
        real, dimension(NXB,NYB,NZB) :: rhoxr,rhoxl,rhoyr,rhoyl,rhozl,rhozr,&
                                        aixr,aixl,aiyr,aiyl,aizr,aizl

        rhoxr = (smrh(ix1:ix2,jy1:jy2,kz1:kz2)+smrh(ix1+1:ix2+1,jy1:jy2,kz1:kz2))/2.0d0
        rhoxl = (smrh(ix1:ix2,jy1:jy2,kz1:kz2)+smrh(ix1-1:ix2-1,jy1:jy2,kz1:kz2))/2.0d0
        rhoyr = (smrh(ix1:ix2,jy1:jy2,kz1:kz2)+smrh(ix1:ix2,jy1+1:jy2+1,kz1:kz2))/2.0d0
        rhoyl = (smrh(ix1:ix2,jy1:jy2,kz1:kz2)+smrh(ix1:ix2,jy1-1:jy2-1,kz1:kz2))/2.0d0

#if NDIM == 3
        rhozr = (smrh(ix1:ix2,jy1:jy2,kz1:kz2)+smrh(ix1:ix2,jy1:jy2,kz1+1:kz2+1))/2.0d0
        rhozl = (smrh(ix1:ix2,jy1:jy2,kz1:kz2)+smrh(ix1:ix2,jy1:jy2,kz1-1:kz2-1))/2.0d0
#endif

        aixr = mdot(ix1:ix2,jy1:jy2,kz1:kz2)*xnorm(ix1:ix2,jy1:jy2,kz1:kz2)*(rhoxr - smrh(ix1:ix2,jy1:jy2,kz1:kz2))
        aixl = mdot(ix1:ix2,jy1:jy2,kz1:kz2)*xnorm(ix1:ix2,jy1:jy2,kz1:kz2)*(rhoxl - smrh(ix1:ix2,jy1:jy2,kz1:kz2))
        aiyr = mdot(ix1:ix2,jy1:jy2,kz1:kz2)*ynorm(ix1:ix2,jy1:jy2,kz1:kz2)*(rhoyr - smrh(ix1:ix2,jy1:jy2,kz1:kz2))
        aiyl = mdot(ix1:ix2,jy1:jy2,kz1:kz2)*ynorm(ix1:ix2,jy1:jy2,kz1:kz2)*(rhoyl - smrh(ix1:ix2,jy1:jy2,kz1:kz2))

#if NDIM == 3
        aizr = mdot(ix1:ix2,jy1:jy2,kz1:kz2)*znorm(ix1:ix2,jy1:jy2,kz1:kz2)*(rhozr - smrh(ix1:ix2,jy1:jy2,kz1:kz2))
        aizl = mdot(ix1:ix2,jy1:jy2,kz1:kz2)*znorm(ix1:ix2,jy1:jy2,kz1:kz2)*(rhozl - smrh(ix1:ix2,jy1:jy2,kz1:kz2))
#endif

        do k=kz1,kz2
                do j=jy1,jy2
                        do i=ix1,ix2        
                                uconv = uconv + mdot(i,j,k)*smrh(i,j,k)
                        end do
                end do
        end do

end subroutine Heat_getConvel
