subroutine Plasma_spReactions(RSP0, RSP1, RSP2, RSP3, RSP4, RSP5, RSP6, RSP7,   &
                              RSP8, RSP9, RSP10, RSP11, RSP12, RSP13, T_h, T_e, & 
                              ix1, ix2, jy1, jy2) 

   ! temperatures have to be in  eV in equations below
   use Plasma_data, ONLY : pls_KtoeV

   implicit none

   real, dimension(:,:,:), intent(inout) :: RSP0,RSP1,RSP2,RSP3,RSP4,&
                                            RSP5,RSP6,RSP7,RSP8,RSP9,&
                                            RSP10,RSP11,RSP12,RSP13
   
   real, dimension(:,:,:), intent(in) :: T_h, T_e
   
   integer, intent(in) :: ix1, ix2, jy1, jy2
   integer :: i,j

   do j=jy1,jy2
     do i=ix1,ix2
        RSP0(i,j,1)  = (1e-6)*(1.5e-9)*((pls_KtoeV*T_e(i,j,1))**0.68)*&
                       EXP(-24.6/(pls_KtoeV*T_e(i,j,1)))
        RSP1(i,j,1)  = (1e-6)*(6e-20)*((T_e(i,j,1)/T_h(i,j,1))**(-4.0))
        RSP2(i,j,1)  = (1e-6)*(5e-10)
        RSP3(i,j,1)  = (1e-6)*(1e-10)*((1.5*pls_KtoeV*T_e(i,j,1))**1.9)*&
                       EXP(-14.6/(1.5*pls_KtoeV*T_e(i,j,1)))
        RSP4(i,j,1)  = (1e-6)*(1.66*1e-6)/((T_h(i,j,1))**0.7)
        RSP5(i,j,1)  = (1e-6)*(1e-12)
        RSP6(i,j,1)  = (1e-6)*(2.03e-8)*((1.5*pls_KtoeV*T_e(i,j,1))**(-0.1))*&
                       EXP(-8.47/(1.5*pls_KtoeV*T_e(i,j,1)))
        RSP7(i,j,1)  = (1e-6)*(2.63e-10)*((1.5*pls_KtoeV*T_e(i,j,1))**(-0.495))*&
                       EXP(-5.65/(1.5*pls_KtoeV*T_e(i,j,1)))
        RSP8(i,j,1)  = (1e-6)*(9.54e-6)*((1.5*pls_KtoeV*T_e(i,j,1))**(-1.05))*&
                       EXP(-55.6/(1.5*pls_KtoeV*T_e(i,j,1)))
        RSP9(i,j,1)  = (1e-6)*(1.68e-5)/(T_h(i,j,1)**0.7)
        RSP10(i,j,1) = (1e-6)*(1.4e-10)
        RSP11(i,j,1) = (1e-6)*(2.6e-10)
        RSP12(i,j,1) = (1e-6)*(1.05e-12)*((pls_KtoeV*T_h(i,j,1))**0.5)
        RSP13(i,j,1) = (1e-6)*(4.5e-12)*&
                       EXP(-3320./(T_h(i,j,1)))
     end do
   end do

end subroutine Plasma_spReactions
