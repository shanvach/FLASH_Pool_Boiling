subroutine Plasma_spGeneration(N_h0,N_h1,N_h2,N_h3,N_h4,N_h5,N_h6,& 
                               N_h7,N_h8,N_h9,N_e,RSP0,RSP1,RSP2, &
                               RSP3,RSP4,RSP5,RSP6,RSP7,RSP8,RSP9,&
                               RSP10,RSP11,RSP12,RSP13,GNH0,GNH1, &
                               GNH2,GNH3,GNH4,GNH5,GNH6,GNH7,GNH8,&
                               GNH9,GNE,GNEBZ,GNERT,ix1,ix2,jy1,jy2)

   implicit none

   real, dimension(:,:,:), intent(in) :: N_h0,N_h1,N_h2,N_h3,N_h4,&
                                         N_h5,N_h6,N_h7,N_h8,N_h9,&
                                         N_e
   
   real, dimension(:,:,:), intent(in) :: RSP0,RSP1,RSP2,RSP3,RSP4,&
                                         RSP5,RSP6,RSP7,RSP8,RSP9,&
                                         RSP10,RSP11,RSP12,RSP13

   real, dimension(:,:,:), intent(inout) :: GNH0,GNH1,GNH2,GNH3,&
                                            GNH4,GNH5,GNH6,GNH7,&
                                            GNH8,GNH9,GNE,GNEBZ,&
                                            GNERT 
   integer, intent(in) :: ix1, ix2, jy1, jy2
   integer :: i,j
   
   do j=jy1,jy2
     do i=ix1,ix2
        !heavy species generation (m3/s)
        GNH0(i,j,1)  = RSP1(i,j,1)*N_e(i,j,1)*N_h1(i,j,1)   +  & 
                       RSP2(i,j,1)*N_h1(i,j,1)*N_h2(i,j,1)  -  &
                       RSP0(i,j,1)*N_e(i,j,1)*N_h0(i,j,1)
        GNH1(i,j,1)  = RSP0(i,j,1)*N_e(i,j,1)*N_h0(i,j,1)   -  &
                       RSP1(i,j,1)*N_e(i,j,1)*N_h1(i,j,1)   -  &
                       RSP2(i,j,1)*N_h1(i,j,1)*N_h2(i,j,1)
        GNH2(i,j,1)  = RSP5(i,j,1)*N_h4(i,j,1)*N_h3(i,j,1)  +  &
                       RSP12(i,j,1)*N_h4(i,j,1)*N_h9(i,j,1) -  &
                       RSP2(i,j,1)*N_h3(i,j,1)*N_h0(i,j,1)  -  &
                       RSP3(i,j,1)*N_e(i,j,1)*N_h2(i,j,1)
        GNH3(i,j,1)  = RSP2(i,j,1)*N_h1(i,j,1)*N_h2(i,j,1)  +  &
                       RSP3(i,j,1)*N_e(i,j,1)*N_h2(i,j,1)   -  &
                       RSP4(i,j,1)*N_e(i,j,1)*N_h3(i,j,1)   -  &
                       RSP5(i,j,1)*N_h3(i,j,1)*N_h4(i,j,1)
        GNH4(i,j,1)  = RSP4(i,j,1)*N_e(i,j,1)*N_h3(i,j,1)   -  &
                       RSP11(i,j,1)*N_h4(i,j,1)*N_h6(i,j,1) -  &
                       RSP12(i,j,1)*N_h4(i,j,1)*N_h9(i,j,1)
        GNH5(i,j,1)  = RSP10(i,j,1)*N_h6(i,j,1)*N_h8(i,j,1) -  &
                       (RSP6(i,j,1)+RSP7(i,j,1)+RSP8(i,j,1))*  &
                       N_e(i,j,1)*N_h5(i,j,1)
        GNH6(i,j,1)  = RSP7(i,j,1)*N_e(i,j,1)*N_h5(i,j,1)   -  &
                       RSP10(i,j,1)*N_h6(i,j,1)*N_h8(i,j,1) -  &
                       RSP11(i,j,1)*N_h4(i,j,1)*N_h6(i,j,1)
        GNH7(i,j,1)  = RSP8(i,j,1)*N_e(i,j,1)*N_h5(i,j,1)   -  &
                       RSP9(i,j,1)*N_e(i,j,1)*N_h7(i,j,1)
        GNH8(i,j,1)  = RSP6(i,j,1)*N_e(i,j,1)*N_h5(i,j,1)   +  &
                       RSP9(i,j,1)*N_e(i,j,1)*N_h7(i,j,1)  -  &
                       RSP10(i,j,1)*N_h6(i,j,1)*N_h8(i,j,1)
        GNH9(i,j,1)  = RSP11(i,j,1)*N_h4(i,j,1)*N_h6(i,j,1) +  &
                       RSP13(i,j,1)*N_h4(i,j,1)*N_h5(i,j,1) -  &
                       RSP12(i,j,1)*N_h4(i,j,1)*N_h9(i,j,1)
        !electrons generation (m3/s)
        GNE(i,j,1)   = RSP0(i,j,1)*N_e(i,j,1)*N_h0(i,j,1)   +  &
                       RSP3(i,j,1)*N_e(i,j,1)*N_h2(i,j,1)   +  &
                       RSP6(i,j,1)*N_e(i,j,1)*N_h5(i,j,1)   +  &
                       RSP8(i,j,1)*N_e(i,j,1)*N_h5(i,j,1)   +  & 
                       RSP11(i,j,1)*N_h4(i,j,1)*N_h6(i,j,1) -  &
                       RSP1(i,j,1)*N_e(i,j,1)*N_h1(i,j,1)    -  &
                       RSP4(i,j,1)*N_e(i,j,1)*N_h2(i,j,1)   -  &
                       RSP7(i,j,1)*N_e(i,j,1)*N_h5(i,j,1)   -  &
                       RSP9(i,j,1)*N_e(i,j,1)*N_h9(i,j,1) 
        !compare with Boltzmann Equilibrium
        GNEBZ(i,j,1) = GNH1(i,j,1) + GNH3(i,j,1) + &
                       GNH7(i,j,1) - GNH6(i,j,1)
        GNERT(i,j,1) = GNEBZ(i,j,1)/GNE(i,j,1)
     end do
   end do 

end subroutine Plasma_spGeneration
