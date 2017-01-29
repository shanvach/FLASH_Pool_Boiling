!!****if* source/Simulation/SimulationMain/DiffuseCtC/Simulation_initSpecies
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
!!  the TwoGamma setup, which advects two fluids with different Gamma values
!!
!!***

subroutine Simulation_initSpecies()
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Multispecies_interface, ONLY : Multispecies_setProperty
!!$  use Simulation_data, ONLY : sim_abar, sim_zbar

  implicit none
#include "Flash.h"
#include "Multispecies.h"

  call RuntimeParameters_get('sim_abar', sim_abar)
  call RuntimeParameters_get('sim_zbar', sim_zbar)

  call Multispecies_setProperty(HE_SPEC, A, sim_abar)
  call Multispecies_setProperty(HE_SPEC, Z, sim_zbar)
  call Multispecies_setProperty(HE_SPEC, GAMMA, 1.6666666666666666667e0)

end subroutine Simulation_initSpecies

