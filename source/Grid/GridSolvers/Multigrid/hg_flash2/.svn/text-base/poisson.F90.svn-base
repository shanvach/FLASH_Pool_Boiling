!*******************************************************************************

!  Routine:     poisson()

!  Description: Driver routine for the Huang-Greengard multigrid Poisson solver.
!               This routine interprets the supplied boundary conditions and
!               either calls the Poisson solver directly (in the case of
!               periodic or Dirichlet boundaries) or uses James' image-mass
!               method to handle isolated boundaries.

!  Parameters:  isoln           Index for solution array.  The solution is
!                               written directly into this variable.
!               isrc            Index for source array.  This is taken to be
!                               the density field; the source array to be used
!                               as the right-hand side of the Poisson equation
!                               is computed from this.
!               bc_types(6)     Boundary condition types array: 
!                                 values = MG_BND_PERIODIC,
!                                          MG_BND_DIRICHLET, or
!                                          MG_BND_ISOLATED 
!                                 index 1 = -x, 2 = +x, 3 = -y, 4 = +y,
!                                       5 = -z, 6 = +z
!               bc_values(2,6)  Values for dirichlet and neumann boundary
!                               conditions.  Not used; present only for
!                               interface compatibility.
!               poisfact        Constant Poisson factor.  Used to scale the
!                               source function.  For example, for gravity,
!                               poisfact = 4*pi*G.


subroutine poisson (isoln, isrc, bc_types, bc_values, poisfact)

!===============================================================================

use mg_common
use dBase, ONLY: dBaseKeyNumber
use perfmon

implicit none

integer               :: isoln, isrc
integer, dimension(6) :: bc_types
real, dimension(2,6)  :: bc_values
real                  :: poisfact

integer               :: i

integer, save         :: isls, icorr, iimg_pot, iimg_mass
logical, save         :: first_call = .true.

external poisson_hg_solve_block

!===============================================================================

call timer_start("poisson")

! Get key numbers from the database for the temporary variables we need.

if (first_call) then
  isls      = dBaseKeyNumber("hgw1")
  icorr     = dBaseKeyNumber("hgw2")
  iimg_mass = dBaseKeyNumber("hgw3")
  iimg_pot  = dBaseKeyNumber("hgw4")
  first_call = .false.
endif

! Perform some generic initializations.

call hg_init()

!-------------------------------------------------------------------------------

! Call the Huang-Greengard Poisson solver for different types of boundary
! conditions.

select case (bc_types(1))

!-------------------------------------------------------------------------------

! In order for isolated boundary conditions to work, an appropriate submodule
! that defines the appropriate mesh variables and implements poisson_image_*
! must be included in the code.

  case (MG_BND_ISOLATED)

    if ((iimg_mass < 1) .or. (iimg_pot < 1)) then
      call abort_flash("poisson:  must include isolated boundary submodule")
    endif

    call hg_solve(isrc, isoln, isls, icorr, poisson_hg_solve_block, &
                  MG_BND_DIRICHLET, poisfact)

    call poisson_image_mass(isoln, iimg_mass)

    do i = 1, mesh_lrefmax
      call hg_level_zero(i, iimg_pot, MG_NODES_ALL_NODES)
    enddo

    call poisson_image_boundary(iimg_mass, iimg_pot, 1.)

    call hg_solve(iimg_mass, iimg_pot, isls, icorr, poisson_hg_solve_block, &
                  MG_BND_GIVENVAL, -1.)

    do i = 1, mesh_lrefmax
      call hg_level_add(i, isoln, iimg_pot, MG_NODES_LEAF_ONLY)
    enddo

!-------------------------------------------------------------------------------

  case (MG_BND_PERIODIC)

    call hg_solve(isrc, isoln, isls, icorr, poisson_hg_solve_block, &
                  MG_BND_PERIODIC, poisfact)

!-------------------------------------------------------------------------------

  case (MG_BND_DIRICHLET)

    call hg_solve(isrc, isoln, isls, icorr, poisson_hg_solve_block, &
                  MG_BND_DIRICHLET, poisfact)

!-------------------------------------------------------------------------------

  case (MG_BND_GIVENVAL)

    call hg_solve(isrc, isoln, isls, icorr, poisson_hg_solve_block, &
                  MG_BND_GIVENVAL, poisfact)

!-------------------------------------------------------------------------------

  case default

    call abort_flash("poisson:  invalid boundary condition type!")

!-------------------------------------------------------------------------------

end select

call timer_stop("poisson")

!===============================================================================

return
end subroutine poisson
