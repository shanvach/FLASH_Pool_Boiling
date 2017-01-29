!!****if* source/Simulation/SimulationMain/unitTest/particleTests/VParticles/Simulation_init
!!
!! NAME
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!  Simulation_init( )
!!
!! DESCRIPTION   
!!   Initialize all the runtime parameters needed for testing virtual 
!!   particles functionality
!!
!! ARGUMENTS
!!
!!   
!!
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none

#include "constants.h"
#include "Flash.h"
  
  
  integer :: i, j, status
 


  call RuntimeParameters_get("sim_initPosX",   sim_initPos(IAXIS))
  call RuntimeParameters_get("sim_initPosY",   sim_initPos(JAXIS))
  call RuntimeParameters_get("sim_initPosZ",   sim_initPos(KAXIS))

  call RuntimeParameters_get("sim_deltaMoveX",   sim_deltaMove(IAXIS))
  call RuntimeParameters_get("sim_deltaMoveY",   sim_deltaMove(JAXIS))
  call RuntimeParameters_get("sim_deltaMoveZ",   sim_deltaMove(KAXIS))

  call RuntimeParameters_get("sim_p_amb",   sim_p_amb )
  call RuntimeParameters_get("sim_rho_amb", sim_rho_amb )
  call RuntimeParameters_get("sim_vx_amb",  sim_vx_amb )
  call RuntimeParameters_get("sim_vx_multiplier",  sim_vx_multiplier )

!! velocity perturbations
  call RuntimeParameters_get("sim_seed",   sim_seed )
  call RuntimeParameters_get("sim_vx_pert",   sim_vx_pert )
  call RuntimeParameters_get("sim_vy_pert",   sim_vy_pert )
  call RuntimeParameters_get("sim_vz_pert",   sim_vz_pert )

  return
end subroutine Simulation_init
