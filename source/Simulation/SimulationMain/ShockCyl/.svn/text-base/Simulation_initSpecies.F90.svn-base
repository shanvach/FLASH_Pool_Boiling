!!****if* source/Simulation/SimulationMain/ShockCyl/Simulation_initSpecies
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
!!  the ShockCyl setup. The two species in use are "Air" and "SF6".
!!
!!***

subroutine Simulation_initSpecies()
  use Multispecies_interface, ONLY : Multispecies_setProperty

implicit none
#include "Multispecies.h"
#include "Flash.h"


  call Multispecies_setProperty(SF6_SPEC, A, 146.)
  call Multispecies_setProperty(SF6_SPEC, Z, 70.)
  call Multispecies_setProperty(SF6_SPEC, GAMMA, 1.09)
! Doesn't seem to matter at initial conditions whether these are set to zero or undefined... that's good!
!  call Multispecies_setProperty(SF6_SPEC, E, 0.)
!  call Multispecies_setProperty(SF6_SPEC, EB, 0.)
!  call Multispecies_setProperty(SF6_SPEC, N, 0.)

  call Multispecies_setProperty(AIR_SPEC, A, 28.66)
  call Multispecies_setProperty(AIR_SPEC, Z, 14.)
  call Multispecies_setProperty(AIR_SPEC, GAMMA, 1.4)
!  call Multispecies_setProperty(AIR_SPEC, E, 0.)
!  call Multispecies_setProperty(AIR_SPEC, EB, 0.)
!  call Multispecies_setProperty(AIR_SPEC, N, 0.)

end subroutine Simulation_initSpecies
