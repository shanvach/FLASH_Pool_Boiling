subroutine Plasma_ColFreq(vei, vea, N_e, T_e, ix1, ix2, jy1, jy2)

   use Plasma_data, ONLY : pls_Ckb, pls_Cme, pls_Cpi, pls_Ce, pls_gam
   implicit none
 
   real, dimension(:,:,:), intent(in) :: N_e, T_e 
   real, dimension(:,:,:), intent(inout) :: vea, vei

   real :: Kel, Lna
   
   integer, intent(in) :: ix1, ix2, jy1, jy2
   integer :: i,j

   do j=jy1,jy2
     do i=ix1,ix2

     Kel  = ((4.0*pls_Cpi*N_e(i,j,1)*pls_Ce**2)/(pls_Ckb*T_e(i,j,1)))**0.5
     Lna = ALOG(4*pls_Ckb*T_e(i,j,1)/(Kel*(pls_Ce**2)*(pls_gam**2))) -      &
           2*ALOG((2.0)**0.5)
     vei(i,j,1) = Lna*(4.0/3.0)*((2*pls_Cpi)**0.5)*N_e(i,j,1)*              &
                  ((pls_Ckb*T_e(i,j,1)/pls_Cme)**0.5)*                      &
                  (pls_Ce**2/(pls_Ckb*T_e(i,j,1)))**2
     vea(i,j,1) = (4.0/3.0)*1e-20*(N_e(i,j,1))*                             &
                  ((8.0*pls_Ckb*T_e(i,j,1)) /(pls_Cpi*pls_Cme))**0.5
     !set limit on collision frequencies
     if (vei(i,j,1).le.1e7) then
      vei(i,j,1) = 1e7
     elseif(vei(i,j,1).ge.1e13) then
      vei(i,j,1) = 1e13
     end if
     !
     if (vea(i,j,1).le.1e7) then
      vea(i,j,1) = 1e7
     elseif(vea(i,j,1).ge.1e13) then
      vea(i,j,1) = 1e13
     end if
     end do
   end do

end subroutine Plasma_ColFreq
