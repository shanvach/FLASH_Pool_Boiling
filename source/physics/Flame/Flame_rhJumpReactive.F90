!!****f* source/physics/Flame/Flame_rhJumpReactive
!!
!! NAME
!!
!!  Flame_rhJumpReactive
!!
!! SYNOPSIS
!!
!!  Flame_rhJumpReactive ( real(:), intent(inout) :: eosData_u,
!!                         real, intent(in) :: qbar_u,
!!                         real(:), intent(out) :: eosData_b,
!!                         real, intent(out) :: qbar_b,
!!                         integer, intent(in) :: eos_mode )
!!
!! DESCRIPTION
!!
!!   Makes a rh reactive jump.
!!
!! ARGUMENTS
!!
!!   eosData_u - equation of state data at u.
!!      qbar_u - ? at u.
!!   eosData_b - equation of state data at b.
!!      qbar_b - ? at b.
!!    eos_mode - equation of state mode.
!!
!!***

subroutine Flame_rhJumpReactive(eosData_u, qbar_u, eosData_b, qbar_b, eos_mode)
  implicit none

#include "Eos.h"

  real, dimension(EOS_NUM), intent(inout) :: eosData_u
  real,    intent(in)                     :: qbar_u
  real, dimension(EOS_NUM), intent(out)   :: eosData_b
  real,    intent(out)                    :: qbar_b
  integer, intent(in)                     :: eos_mode

  eosData_b = 0.0
  qbar_b = 0.0
end subroutine Flame_rhJumpReactive
