!*******************************************************************************

! Routine:      poisson_image_boundary

! Description:  Compute the boundary values of the image mass potential.


subroutine poisson_image_boundary (iiden, iipot, poisfact)

!===============================================================================

implicit none

integer       :: iiden, iipot
real          :: poisfact
logical, save :: first_call = .true.

!===============================================================================

if (first_call) then
  call init_mpole
  first_call = .false.
endif

call find_center_of_mass (iiden)
call compute_mpole_moments (iiden)
call compute_mpole_potential (iipot, poisfact)

!===============================================================================

return
end
