!!****if* source/Simulation/SimulationMain/RayTraceGradient/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Initializes initial conditions for Sedov Spherical Explosion 
!!  problem.
!!
!! ARGUMENTS
!!
!!   
!!
!! PARAMETERS
!!
!!  p_ambient       Initial ambient pressure
!!  rho_Ambient     Initial ambient density
!!  exp_energy      Explosion energy (distributed over 2^dimen central zones)
!!  t_init          Initial time since explosion
!!  sim_nsubzones      Number of `sub-zones' in cells for applying 1d profile
!!
!!***
subroutine Simulation_init()

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "constants.h"
#include "Flash.h"

  call RuntimeParameters_get('sim_temp', sim_temp)
  call RuntimeParameters_get('sim_gradient', sim_gradient)
  call RuntimeParameters_get('sim_ncrit', sim_ncrit)

  return
end subroutine Simulation_init
