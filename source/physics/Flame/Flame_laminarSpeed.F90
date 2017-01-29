!!****f* source/physics/Flame/Flame_laminarSpeed
!!
!! NAME
!!
!!  Flame_laminarSpeed
!!
!! SYNOPSIS
!!
!!  call Flame_laminarSpeed (real,intent(in) :: dens,
!!                           real,intent(in) :: pres,
!!                           real,intent(out) :: s )
!!
!! DESCRIPTION
!!
!!   Determines the laminar speed.
!!
!! ARGUMENTS
!!
!!   dens - density
!!   pres - pressure
!!   s    - speed
!!
!!***

!------------------------------------------------------------------------

subroutine Flame_laminarSpeed(dens, pres, s)


  implicit none

  real, intent(in)   :: dens
  real, intent(in)   :: pres
  real, intent(out)  :: s

  s = 0.0
  
end subroutine Flame_laminarSpeed

!------------------------------------------------------------------------

