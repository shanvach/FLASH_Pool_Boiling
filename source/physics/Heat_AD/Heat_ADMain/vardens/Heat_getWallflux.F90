subroutine Heat_getWallflux(pf,T,Nu_l,Nu_t,hcounter,dy,ycell,jy1,ix1,ix2,kz1,kz2)

     implicit none
     real, dimension(:,:,:),intent(in) :: pf, T
     real, intent(inout) :: Nu_l,Nu_t
     integer, intent(inout) :: hcounter
     real, intent(in) :: dy,ycell
     integer, intent(in) :: jy1,ix1,ix2,kz1,kz2

     integer :: i,k

     if(ycell == 0.5*dy) then

      do k=kz1,kz2
       do i=ix1,ix2

          Nu_l = Nu_l + (1.0 - pf(i,jy1,k))*(1.0 - T(i,jy1,k))/(0.5*dy)
          Nu_t = Nu_t + (1.0 -  T(i,jy1,k))/(0.5*dy)
          hcounter = hcounter + 1

        end do
      end do             

     end if

end subroutine Heat_getWallflux
