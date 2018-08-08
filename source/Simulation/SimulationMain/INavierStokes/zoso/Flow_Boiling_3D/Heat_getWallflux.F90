subroutine Heat_getWallflux(pf,T,Hf,Nu_l,Nu_t,hcounter,dy,ycell,jy1,ix1,ix2,kz1,kz2,blockID)

     use Grid_interface, ONLY : Grid_getBlkBoundBox, Grid_getBlkCenterCoords, Grid_getDeltas
     use Multiphase_data,ONLY : mph_thco1,mph_thco2

#include "Flash.h"
#include "constants.h"

     implicit none
     real, dimension(:,:,:),intent(in)    :: pf, T
     real, dimension(:,:,:),intent(inout) :: Hf
     integer, intent(in) :: blockID
     real, intent(inout) :: Nu_l,Nu_t
     integer, intent(inout) :: hcounter
     real, intent(in) :: dy,ycell
     integer, intent(in) :: jy1,ix1,ix2,kz1,kz2

     integer :: i,k
     real del(MDIM),bsize(MDIM),coord(MDIM)
     real, dimension(2,MDIM) :: boundBox
     real :: xcell, zcell
     real :: tol=1E-13

     call Grid_getDeltas(blockID,del)
     call Grid_getBlkCenterCoords(blockId,coord)
     call Grid_getBlkBoundBox(blockId,boundBox)

     bsize(:) = boundBox(2,:) - boundBox(1,:)

     if(abs(ycell-0.5*dy) .le. tol) then

      do k=kz1,kz2
       do i=ix1,ix2

          xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
                  real(i - NGUARD - 1)*del(IAXIS) +   &
                  0.5*del(IAXIS)

          zcell  = coord(KAXIS) - bsize(KAXIS)/2.0 +  &
                   real(k - NGUARD - 1)*del(KAXIS)  +  &
                   0.5*del(KAXIS)

          if(xcell .ge. -4.95 .and. xcell .le. 4.95) then

                Nu_l = Nu_l + (1.0 - pf(i,jy1,k))*(1.0 - T(i,jy1,k))/(0.5*dy)
                Nu_t = Nu_t + ((1.0 - T(i,jy1,k))/(0.5*dy))*(1.0-pf(i,jy1,k)+(mph_thco1/mph_thco2)*pf(i,jy1,k))
                hcounter = hcounter + 1

                Hf(i,jy1,k) = (1.0 - pf(i,jy1,k))*(1.0 - T(i,jy1,k))/(0.5*dy)

          end if

        end do
      end do             

     end if

end subroutine Heat_getWallflux
