subroutine Heat_Solve(T, Trhs, Told, dt, ix1, ix2, jy1, jy2, kz1, kz2, gama, rho)

  implicit none

  real, dimension(:,:,:), intent(inout) :: T
  real, dimension(:,:,:), intent(in)    :: Trhs, Told
  real, intent(in) :: dt, gama, rho
  integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2

  T(ix1:ix2,jy1:jy2,kz1:kz2) =                T(ix1:ix2,jy1:jy2,kz1:kz2) &
                             + dt * gama * Trhs(ix1:ix2,jy1:jy2,kz1:kz2) &
                             + dt * rho  * Told(ix1:ix2,jy1:jy2,kz1:kz2)

end subroutine Heat_Solve
