!!****f* source/physics/Flame/Flame_rhJump
!!
!! NAME
!!
!!  Flame_rhJump
!!
!! SYNOPSIS
!!
!!  call Flame_rhJump (
!!    real, intent(inout) :: dens_u,
!!    real, intent(inout) :: pres_u,
!!    real, intent(inout) :: temp_u,
!!    real, intent(inout) :: ener_u,
!!    real, intent(in) :: ye_u,
!!    real, intent(in) :: sumy_u,
!!    real, intent(out) :: dens_b,
!!    real, intent(out) :: pres_b,
!!    real, intent(out) :: temp_b,
!!    real, intent(out) :: ener_b,
!!    real, intent(in) :: ye_b,
!!    real, intent(in) :: sumy_b,
!!    real, intent(in) :: q,
!!    real, intent(in) :: s,
!!    integer, intent(in) :: flag )
!!
!! DESCRIPTION
!!
!!   Makes a rh jump.
!!
!! ARGUMENTS
!!
!!   dens_u - density at u.
!!   pres_u - pressure at u.
!!   temp_u - temperature at u.
!!   ener_u - energy at u.
!!   ye_u   - ? at u.
!!   sumy_u - ? at u.
!!   dens_b - density at b.
!!   pres_b - pressure at b.
!!   temp_b - temperature at b.
!!   ener_b - energy at b.
!!   ye_b   - ? at b.
!!   sumy_b - ? at b.
!!   q      - ?
!!   s      - ?
!!   flag   - ?
!!
!!***

subroutine Flame_rhJump(dens_u, pres_u, temp_u, ener_u, ye_u, sumy_u, &
          dens_b, pres_b, temp_b, ener_b, ye_b, sumy_b, q, s, flag)

  implicit none
       
  real,    intent(inout)  :: dens_u, pres_u, temp_u, ener_u
  real,    intent(out)    :: dens_b, pres_b, temp_b, ener_b
  real,    intent(in)     :: q, s
  integer, intent(in)     :: flag
  real, intent(in) ::  ye_u, sumy_u, ye_b, sumy_b

  dens_b=dens_u
  pres_b=pres_u
  temp_b=temp_u
  ener_b=pres_u
  return
end subroutine Flame_rhJump
