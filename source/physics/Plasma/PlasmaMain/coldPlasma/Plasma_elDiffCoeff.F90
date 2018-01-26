subroutine Plasma_elDiffCoeff(Diffcoeffion, DiffCoeffel, T_e, T_h,& 
                              Fvea, Fvei, ix1, ix2, jy1, jy2)
   
   use Plasma_data, ONLY : pls_Ckb, pls_Cmi_net

   real, dimension(:,:,:), intent(inout) :: DiffCoeffion,DiffCoeffel
   real, dimension(:,:,:), intent(in) :: T_e, T_h, Fvea, Fvei

   integer, intent(in) :: ix1, ix2, jy1, jy2
   integer :: i,j

   do j=jy1,jy2
     do i=ix1,ix2
     !net ion coefficient
     DiffCoeffion(i,j,1) = ( pls_Ckb*T_h(i,j,1))/&
                           ( pls_Cmi_net*(Fvea(i,j,1)+Fvei(i,j,1)))
     !effective ion and electron coefficient
     Diffcoeffel(i,j,1)  = DiffCoeffion(i,j,1)*(1.0 + T_e(i,j,1)/T_h(i,j,1))
     end do
   end do

end subroutine Plasma_elDiffCoeff
