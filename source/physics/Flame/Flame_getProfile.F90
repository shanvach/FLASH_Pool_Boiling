!!****f* source/physics/Flame/Flame_getProfile
!!
!! NAME
!!
!!  Flame_getProfile
!!
!! SYNOPSIS
!!
!!  call Flame_getProfile (real,intent(in) :: x, real,intent(out) :: f)
!!
!! DESCRIPTION
!!   Gets the profile
!!
!! ARGUMENTS
!!
!!   x - ?
!!   f - ?
!!
!!***

subroutine Flame_getProfile(x, f)

  implicit none

  real, intent(in)  :: x
  real, intent(out) :: f

  f = 0.0

end subroutine Flame_getProfile
