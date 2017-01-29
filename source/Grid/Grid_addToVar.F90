!!****f* source/Grid/Grid_addToVar
!!
!! NAME
!!
!!  Grid_addToVar
!!
!! SYNOPSIS
!!
!!  call Grid_addToVar(integer(in) :: srcVar(:),
!!                     integer(in) :: destVar(:),
!!                     real(in) :: multFactor(:),
!!                     logical(in) :: reset)
!!
!! DESCRIPTION
!!   Compute solnData(srcVar,:,:,:)*multFactor 
!!   to solnData(destVar,:,:,:)
!!   If reset is true, the target is first zeroed
!!   For a copy call with multFactor=1.0 and reset=.true.
!!   srcVar == destVar is allowed
!!
!!
!! ARGUMENTS
!!
!!
!!   srcVar : the state variables to be used in the RHS of the expression
!!
!!   destVar : the state variables to be used in the LHS of the expression
!!
!!   multFactor : multiplication factor
!!
!!   reset : indicates whether the destination variable should be zeroed first
!!
!!
!!
!!***

subroutine Grid_addToVar(srcVar, destVar, multFactor, reset)

  implicit none
  integer, intent(in) :: srcVar, destVar
  real,  intent(in) :: multFactor
  logical, intent(in) :: reset

end subroutine Grid_addToVar
