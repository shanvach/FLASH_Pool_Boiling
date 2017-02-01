subroutine Heat_calMdot(mdot,Tnl,Tnv,alpha_l,alpha_v,nx,ny,ix1,ix2,jy1,jy2,kz1,kz2)
        
     use Heat_AD_data
     use IncompNS_data, only: ins_invRe

     implicit none
     real, dimension(:,:,:), intent(inout) :: mdot
     real, dimension(:,:,:), intent(in) :: Tnl,Tnv,nx,ny
     integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
     real, intent(in) :: alpha_l, alpha_v

     integer :: i,j,k
     real :: SRP,alpha

     SRP = (ht_St*ins_invRe)/(ht_Pr)
     alpha = alpha_v/alpha_l    

     do k=kz1,kz2
      do j=jy1,jy2
       do i=ix1,ix2

        mdot(i,j,k) = SRP*(Tnl(i,j,k) + alpha*Tnv(i,j,k))
        !mdot(i,j,k) = ((alpha_l*(Tnl(i,j,k)) + alpha_v*(Tnv(i,j,k))))/ht_L
        !mdot(i,j,k) = -0.1
 
       end do
     end do
    end do

end subroutine Heat_calMdot
