subroutine Heat_Solve2d(T_p, T_o, u, v, dt, dx, dy, ix1,ix2, jy1, jy2)

  implicit none

  real, dimension(:,:,:), intent(inout) :: T_p
  real, dimension(:,:,:), intent(in) :: T_o
  real, dimension(:,:,:), intent(in) :: u,v
  real, intent(in) :: dt
  real, dimension(:,:),   intent(in) :: dx, dy 
  integer, intent(in) :: ix1, ix2, jy1, jy2

end subroutine Heat_Solve2d


subroutine Heat_Solve3d(T_p, T_o, u, v, w, dt, dx, dy, dz, &
                        ix1,ix2, jy1, jy2, kz1, kz2)

  implicit none

  real, dimension(:,:,:), intent(inout) :: T_p
  real, dimension(:,:,:), intent(in) :: T_o
  real, dimension(:,:,:), intent(in) :: u,v,w
  real, intent(in) :: dt
  real, dimension(:,:),   intent(in) :: dx, dy, dz 
  integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2

end subroutine Heat_Solve3d
