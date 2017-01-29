!!****f* source/physics/Flame/Flame_getWidth
!!
!! NAME
!!
!!  Flame_getWidth
!!
!! SYNOPSIS
!!
!!  call Flame_getWidth ( real, intent(out) :: laminarWidth )
!!
!! DESCRIPTION
!!
!!  Gets the laminar width.
!!
!! ARGUMENTS
!!
!!  laminarWidth  -  the laminar width.
!!
!!***

subroutine Flame_getWidth(laminarWidth)

  implicit none

  real, intent(OUT) :: laminarWidth

  laminarWidth = 0.0

end subroutine Flame_getWidth

