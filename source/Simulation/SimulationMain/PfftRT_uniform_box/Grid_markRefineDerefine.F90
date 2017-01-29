!!****if* source/Simulation/SimulationMain/PfftRT_uniform_box/Grid_markRefineDerefine
!!
!! NAME
!!  Grid_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  Grid_markRefineDerefine()
!!  
!! DESCRIPTION 
!!  Mark blocks for refinement or derefinement
!!  This routine is used with AMR only where individual 
!!  blocks are marked for refinement or derefinement based upon
!!  some refinement criterion. The Uniform Grid does not need
!!  this routine, and uses the stub.
!!
!! ARGUMENTS
!! 
!! NOTES
!!
!! Every unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. For Grid unit these variables begin with "gr_"
!! like, gr_meshMe or gr_eosMode, and are stored in fortran
!! module Grid_data (in file Grid_data.F90). The other variables
!! are local to the specific routines and do not have the prefix "gr_"
!!
!!
!!***

subroutine Grid_markRefineDerefine()

  use tree, ONLY : lrefine_max
  use Simulation_data
  use Grid_interface, ONLY : Grid_markRefineSpecialized, Grid_getDomainBoundBox

  implicit none

#include "constants.h"
#include "Flash.h"

  
  real, dimension(2, MDIM) :: region_bb
  real, dimension(7) :: rect_specs
  real :: x_ctr, x_wid
  
  call Grid_getDomainBoundBox(region_bb)
  x_ctr = 0.5 * ( region_bb(2, IAXIS) + region_bb(1, IAXIS) )
  x_wid = 0.5 * ( region_bb(2, IAXIS) - region_bb(1, IAXIS) )
  rect_specs(1) = x_ctr - x_wid * xrefine_rel
  rect_specs(2) = x_ctr + x_wid * xrefine_rel
  rect_specs(3) = region_bb(1, JAXIS)
  rect_specs(4) = region_bb(2, JAXIS)
  rect_specs(5) = region_bb(1, KAXIS)
  rect_specs(6) = region_bb(2, KAXIS)
  rect_specs(7) = 0
  
  call Grid_markRefineSpecialized(RECTANGLE, 7, rect_specs, lrefine_max)
  
  return
end subroutine Grid_markRefineDerefine

