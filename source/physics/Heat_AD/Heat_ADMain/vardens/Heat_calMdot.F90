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

     ! Mass Flux Equation
    ! mdot(ix1:ix2,jy1:jy2,kz1:kz2) = SRP*(Tnl(ix1:ix2,jy1:jy2,kz1:kz2)+alpha*Tnv(ix1:ix2,jy1:jy2,kz1:kz2))

     ! Constant Mass Flux
     mdot(ix1:ix2,jy1:jy2,kz1:kz2) = 0.0

end subroutine Heat_calMdot
