subroutine Plasma_elDiffCoeff(DiffCoeffel, T_e, Fvea, Fvei, ix1, ix2, jy1, jy2)
   
   use Plasma_data, ONLY : pls_Ckb, pls_Cme

   real, dimension(:,:,:), intent(inout) :: DiffCoeffel
   real, dimension(:,:,:), intent(in) :: T_e, Fvea, Fvei

   integer, intent(in) :: ix1, ix2, jy1, jy2
   integer :: i,j

   do j=jy1,jy2
     do i=ix1,ix2

     DiffCoeffel(i,j,1) = ( pls_Ckb*T_e(i,j,1))/&
                          ( pls_Cme*(Fvea(i,j,1)+Fvei(i,j,1)))
   
     end do
   end do

end subroutine Plasma_elDiffCoeff
