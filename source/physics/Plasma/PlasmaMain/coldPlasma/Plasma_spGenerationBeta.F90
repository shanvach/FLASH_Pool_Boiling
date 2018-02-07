subroutine Plasma_spGeneration(case_num,&
                               N_h0,N_h1,N_h2,N_h3,N_h4,N_h5,&
                               N_h6,N_h7,N_h8,N_h9,N_el,GNRS,&
                               ix1,ix2,jy1,jy2)

   implicit none

   real, dimension(:,:,:), intent(in) :: N_h0,N_h1,N_h2,N_h3,N_h4,&
                                         N_h5,N_h6,N_h7,N_h8,N_h9,&
                                         N_el
  
   real, dimension(:,:,:), intent(inout) :: GNRS
   
   integer, intent(in) :: ix1, ix2, jy1, jy2
   integer :: i,j, case_num

   select case(case_num)

      case(0) !He

         GNRS(i,j,1)  = N_el(i,j,1)*N_h6(i,j,1)*((1e-6)*(6e-20)*((T_e(i,j,1)/T_h(i,j,1))**(-4.0))) +  &
                        N_h1(i,j,1)*N_h6(i,j,1)*((1e-6)*(5e-10))                                   -  &
                        N_el(i,j,1)*N_h0(i,j,1)*((1e-6)*(1.5e-9)*((pls_KtoeV*T_e(i,j,1))**0.68)*&
                                                        EXP(-24.6/(pls_KtoeV*T_e(i,j,1))))
      
      case(1) !N2
      
         GNH1(i,j,1)  = RSP6(i,j,1)*N_h3(i,j,1)*N_h7(i,j,1)   +  &
                        RSP15(i,j,1)*N_h3(i,j,1)*N_h5(i,j,1)  -  &
                        RSP2(i,j,1)*N_h1(i,j,1)*N_h6(i,j,1)   -  &
                        RSP4(i,j,1)*N_el(i,j,1)*N_h1(i,j,1)

      case(2) !O2
    
      case(3) !N 

      case(4) !O
    
      case(5) !NO

      case(6) !He+

      case(7) !N2+

      case(8) !O2+
 
      case(9) !O-

      case(10) !e

   end select
   
   
   do j=jy1,jy2
     do i=ix1,ix2
        !heavy species generation (m3/s)
        GNH0(i,j,1)  = RSP1(i,j,1)*N_el(i,j,1)*N_h6(i,j,1)   +  & 
                       RSP2(i,j,1)*N_h1(i,j,1)*N_h6(i,j,1)   -  &
                       RSP0(i,j,1)*N_el(i,j,1)*N_h0(i,j,1)
        GNH1(i,j,1)  = RSP6(i,j,1)*N_h3(i,j,1)*N_h7(i,j,1)   +  &
                       RSP15(i,j,1)*N_h3(i,j,1)*N_h5(i,j,1)  -  &
                       RSP2(i,j,1)*N_h1(i,j,1)*N_h6(i,j,1)   -  &
                       RSP4(i,j,1)*N_el(i,j,1)*N_h1(i,j,1)
        GNH2(i,j,1)  = RSP11(i,j,1)*N_h4(i,j,1)*N_h9(i,j,1)  -  &
                       ((RSP7(i,j,1)+RSP8(i,j,1)+RSP9(i,j,1))*  &
                         N_el(i,j,1)*N_h2(i,j,1))            -  &
                       RSP16(i,j,1)*N_h2(i,j,1)*N_h3(i,j,1)
        GNH3(i,j,1)  = RSP5(i,j,1)*N_el(i,j,1)*N_h7(i,j,1)   -  &
                       RSP6(i,j,1)*N_h3(i,j,1)*N_h7(i,j,1)   -  &
                       RSP14(i,j,1)*N_h3(i,j,1)*N_h9(i,j,1)  -  &
                       RSP15(i,j,1)*N_h3(i,j,1)*N_h5(i,j,1)  -  &
                       RSP16(i,j,1)*N_h2(i,j,1)*N_h3(i,j,1)
        GNH4(i,j,1)  = RSP7(i,j,1)*N_el(i,j,1)*N_h2(i,j,1)   +  &
                       RSP8(i,j,1)*N_el(i,j,1)*N_h2(i,j,1)   +  &
                       RSP10(i,j,1)*N_el(i,j,1)*N_h8(i,j,1)  +  &
                       RSP16(i,j,1)*N_h2(i,j,1)*N_h3(i,j,1)  -  &
                       RSP11(i,j,1)*N_h4(i,j,1)*N_h9(i,j,1)
        GNH5(i,j,1)  = RSP14(i,j,1)*N_h3(i,j,1)*N_h9(i,j,1)  +  &
                       RSP16(i,j,1)*N_h2(i,j,1)*N_h3(i,j,1)  -  &
                       RSP15(i,j,1)*N_h3(i,j,1)*N_h5(i,j,1) 
        GNH6(i,j,1)  = RSP0(i,j,1)*N_el(i,j,1)*N_h0(i,j,1)   -  &
                       RSP1(i,j,1)*N_el(i,j,1)*N_h6(i,j,1)   -  &
                       RSP2(i,j,1)*N_h1(i,j,1)*N_h6(i,j,1)
        GNH7(i,j,1)  = RSP2(i,j,1)*N_h1(i,j,1)*N_h6(i,j,1)   +  &
                       RSP4(i,j,1)*N_el(i,j,1)*N_h1(i,j,1)   -  &
                       RSP5(i,j,1)*N_el(i,j,1)*N_h7(i,j,1)   -  &
                       RSP6(i,j,1)*N_h3(i,j,1)*N_h7(i,j,1)
        GNH8(i,j,1)  = RSP9(i,j,1)*N_el(i,j,1)*N_h2(i,j,1)   -  &
                       RSP10(i,j,1)*N_el(i,j,1)*N_h8(i,j,1)
        GNH9(i,j,1)  = RSP8(i,j,1)*N_el(i,j,1)*N_h2(i,j,1)   -  &
                       RSP11(i,j,1)*N_h4(i,j,1)*N_h9(i,j,1)  -  &
                       RSP14(i,j,1)*N_h3(i,j,1)*N_h9(i,j,1)
        !electrons generation (m3/s)
        GNE(i,j,1)   = RSP0(i,j,1)*N_el(i,j,1)*N_h0(i,j,1)   +  &
                       RSP4(i,j,1)*N_el(i,j,1)*N_h1(i,j,1)   +  &
                       RSP9(i,j,1)*N_el(i,j,1)*N_h2(i,j,1)   +  &
                       RSP11(i,j,1)*N_h4(i,j,1)*N_h9(i,j,1)  +  & 
                       RSP14(i,j,1)*N_h3(i,j,1)*N_h9(i,j,1)  +  &
                       RSP16(i,j,1)*N_h2(i,j,1)*N_h3(i,j,1)  +  &
                       RSP1(i,j,1)*N_el(i,j,1)*N_h6(i,j,1)   -  &
                       RSP5(i,j,1)*N_el(i,j,1)*N_h7(i,j,1)   -  &
                       RSP8(i,j,1)*N_el(i,j,1)*N_h2(i,j,1)   -  &
                       RSP10(i,j,1)*N_el(i,j,1)*N_h8(i,j,1)  
        !compare with Boltzmann Equilibrium
        GNEBZ(i,j,1) = GNH6(i,j,1) + GNH7(i,j,1) + &
                       GNH8(i,j,1) - GNH9(i,j,1)
        GNERT(i,j,1) = GNEBZ(i,j,1)/GNE(i,j,1)
     end do
   end do 

end subroutine Plasma_spGeneration
