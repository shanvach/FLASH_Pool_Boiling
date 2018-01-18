subroutine Plasma_hvDiffCoeff( DiffCoeff, P_h, T_h, ix1, ix2, jy1, jy2, RSCD, MHSP)
  
   use Plasma_data, ONLY : pls_MMIX
 
   implicit none

   real, dimension(:,:,:), intent(inout) :: DiffCoeff
   real, dimension(:,:,:), intent(in) :: P_h, T_h
   real, intent(in) :: RSCD, MHSP

   integer, intent(in) :: ix1, ix2, jy1, jy2
   integer :: i,j

   do j=jy1,jy2
     do i=ix1,ix2
       
     DiffCoeff(i,j,1) = (1e-4)*((2.63e-7)/(( P_h(i,j,1)/101325.0)*( RSCD )**2))* &
                         (( ((T_h(i,j,1))**3)*( MHSP + pls_MMIX )/ &
                         (2*MHSP*pls_MMIX)))**0.5
     end do
  end do
end subroutine Plasma_hvDiffCoeff

