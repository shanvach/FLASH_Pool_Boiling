!!****if* source/Simulation/SimulationMain/magnetoHD/BubbleShock/Simulation_initSpecies
!!
!! NAME
!!
!!  Simulation_initSpecies
!!
!!
!! SYNOPSIS
!!  Simulation_initSpecies()
!!
!! DESCRIPTION
!!
!!  This routine will initialize the species and species values needed for
!!  the BubbleShock setup, which advects two fluids with different Gamma values
!!
!!***

subroutine Simulation_initSpecies()
  use Multispecies_interface, ONLY : Multispecies_setProperty

  implicit none
#include "Flash.h"
#include "Multispecies.h"

  call Multispecies_setProperty(FLD1_SPEC, A, 1.)
  call Multispecies_setProperty(FLD1_SPEC, Z, 1.)
  call Multispecies_setProperty(FLD1_SPEC, GAMMA, 1.66666666667e0)

  call Multispecies_setProperty(FLD2_SPEC, A, 4.0)
  call Multispecies_setProperty(FLD2_SPEC, Z, 2.0)
  call Multispecies_setProperty(FLD2_SPEC, GAMMA, 2.0)

end subroutine Simulation_initSpecies

