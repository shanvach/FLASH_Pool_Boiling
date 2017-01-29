!*******************************************************************************

! Routine:      poisson_hg_solve_block

! Description:  Single-block Poisson solver for use with the Huang-Greengard
!               multilevel adaptive solver.

!               This version is a stub, meant to be overridden by routines
!               supplied by submodules.

! Parameters:   soln(:,:,:)   Array to receive solution, sampled at zone
!                             centers.  On input contains source.
!               nx, ny, nz    Dimensions of source/solution array
!               dx, dy, dz    Zone widths in x, y, and z
!               bnd_type      Type of boundary conditions to assume
!               level         Level we're solving on


subroutine poisson_hg_solve_block (soln, nx, ny, nz, dx, dy, dz, bnd_type,&
                                   level)

!===============================================================================

implicit none

integer, intent(in) :: nx, ny, nz, bnd_type, level
real, intent(in)    :: dx, dy, dz
real, intent(inout) :: soln(nx,ny,nz)

!===============================================================================

!===============================================================================

return
end subroutine poisson_hg_solve_block
