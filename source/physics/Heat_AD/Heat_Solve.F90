subroutine Heat_Solve(T_p, T_o, u, v, dt, dx, dy, dz, inRe, ix1,ix2, jy1, jy2)

         implicit none
         real, dimension(:,:,:), intent(inout) :: T_p
         real, dimension(:,:,:), intent(in) :: T_o
         real, dimension(:,:,:), intent(in) :: u,v
         real, intent(in) :: dt, dx, dy, dz, inRe
         integer, intent(in) :: ix1, ix2, jy1, jy2

end subroutine Heat_Solve
