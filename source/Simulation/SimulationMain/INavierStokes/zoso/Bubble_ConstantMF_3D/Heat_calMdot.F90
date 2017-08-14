subroutine Heat_calMdot(mdot,Tnl,Tnv,alpha_l,alpha_v,nx,ny,ix1,ix2,jy1,jy2,kz1,kz2)
        
     use Heat_AD_data
     use IncompNS_data, only: ins_invRe

     use RuntimeParameters_interface, ONLY : RuntimeParameters_get

     implicit none
     real, dimension(:,:,:), intent(inout) :: mdot
     real, dimension(:,:,:), intent(in) :: Tnl,Tnv,nx,ny
     integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
     real, intent(in) :: alpha_l, alpha_v

     integer :: i,j,k
     real :: SRP,alpha,beta


     call RuntimeParameters_get("beta",beta)

     ! Constant Mass Flux
     mdot(ix1:ix2,jy1:jy2,kz1:kz2) = -beta

end subroutine Heat_calMdot
