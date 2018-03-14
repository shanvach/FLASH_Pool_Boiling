subroutine Plasma_getNorm(nrmx,nrmy,s,ix1,ix2,jy1,jy2,dx,dy)

     implicit none

     !Arguements
     real, intent(inout), dimension(:,:,:) :: nrmx,nrmy
     real, intent(in), dimension(:,:,:) :: s
     integer, intent(in) :: ix1,ix2,jy1,jy2
     real, intent(in) :: dx,dy

     integer :: kz1 = 1

     nrmx(ix1:ix2,jy1:jy2,kz1) =           &
           -(( s(ix1+1:ix2+1,jy1:jy2,kz1) -   &
             s(ix1-1:ix2-1,jy1:jy2,kz1) )/2./dx)/ &
             sqrt( ((s(ix1+1:ix2+1,jy1:jy2,kz1) - &
             s(ix1-1:ix2-1,jy1:jy2,kz1))/2./dx)**2 &
             + ((s(ix1:ix2,jy1+1:jy2+1,kz1) - &
             s(ix1:ix2,jy1-1:jy2-1,kz1))/2./dy)**2 )

     nrmy(ix1:ix2,jy1:jy2,kz1) =           &
           -(( s(ix1:ix2,jy1+1:jy2+1,kz1) -   &
             s(ix1:ix2,jy1-1:jy2-1,kz1) )/2./dy)/ &
             sqrt( ((s(ix1+1:ix2+1,jy1:jy2,kz1) - &
             s(ix1-1:ix2-1,jy1:jy2,kz1))/2./dx)**2 &
             + ((s(ix1:ix2,jy1+1:jy2+1,kz1) - &
             s(ix1:ix2,jy1-1:jy2-1,kz1))/2./dy)**2 )

end subroutine Plasma_getNorm
