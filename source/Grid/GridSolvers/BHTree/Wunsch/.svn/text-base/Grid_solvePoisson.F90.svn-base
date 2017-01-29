!!****if* source/Grid/GridSolvers/BHTree/Wunsch/Grid_solvePoisson
!!
!! NAME
!!
!!  Grid_solvePoisson
!!
!!
!! SYNOPSIS
!!
!!   Grid_solvePoisson(integer(IN) :: ipotvar,
!!           integer(IN) :: idensvar, 
!!           integer(6)(IN) :: bcTypes,
!!           real(2,6)(IN) :: bcValues,
!!           real(INOUT) :: poisfact)
!!
!! DESCRIPTION
!!
!!   Poisson solver routine.  This module implements the tree
!!   summation method for isolated problems.  Periodic problems are
!!   not supported; 
!!
!!
!! ARGUMENTS
!!
!!  ipotvar -  index to variable containing potential
!!  idensvar - index to variable containing density
!!  bcTypes - boundary types along various faces,
!!             only used in verifying that they are isolated
!!  bcValues - the values to boundary conditions, currently not used
!!  poisfact -  factor to be used in calculation, for gravity: 4*PI*G
!!
!!***



subroutine Grid_solvePoisson (ipotvar, idensvar, bcTypes, bcValues, poisfact)

  use tree, ONLY: grid_changed
  use gr_bhData, ONLY: gr_bhIlist, gr_bhBndType, gr_bhGravFac
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use gr_bhInterface, ONLY : gr_bhComBlkProperties, gr_bhFindNeighbours, &
    gr_bhBuildTree, gr_bhExchangeTrees, gr_bhPotential, gr_bhDestroyTree

  implicit none
#include "constants.h"

  integer, intent(in)    :: ipotvar, idensvar
  integer, intent(in)    :: bcTypes(6)
  real, intent(in)       :: bcValues(2,6)
  real, intent(inout)    :: poisfact

  
  !=======================================================================
  

  call Timers_start("gravity_tree")

  gr_bhBndType = bcTypes
  gr_bhGravFac = poisfact / (4*PI)

  ! if treeBlkProperties was not called yet (start or restart) call it now
  ! later it is called by amr_refine_derefine
  if (grid_changed .eq. 1) then
    call gr_bhComBlkProperties()
    if (gr_bhIlist .eq. 1) call gr_bhFindNeighbours()
  endif

  call gr_bhBuildTree(idensvar)
  call gr_bhExchangeTrees()
  call gr_bhPotential(idensvar, ipotvar)
  call gr_bhDestroyTree()

  call Timers_stop("gravity_tree")
  !=========================================================================
  
  return
end subroutine Grid_solvePoisson
