!!****if* source/Simulation/SimulationMain/unitTest/Particles/Simulation_init
!!
!! NAME
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!  Simulation_init( )
!!
!! DESCRIPTION   
!!   Initialize all the runtime parameters needed for the Particle unitTest
!!
!! ARGUMENTS
!!
!!   
!!
!! PARAMETERS
!!
!!   sim_rho_amb    Gas Density:  Entire domain receives this ambient parameter
!!   sim_p_amb      Gas Pressure: Entire domain receives this ambient parameter
!!   sim_vx_amb     Gas x-velocity:  Dominant flow velocity throughout domain 
!!   sim_vx_multiplier   Half of the domain in y has x-velocity multiplied by this value
!!   sim_seed   Random number seed -- NOT USED please ignore
!!   sim_vx_pert   Scales [-1,1] random number in x direction: set to zero for uniform flow
!!   sim_vy_pert   Scales [-1,1] random number in y direction: set to zero for uniform flow
!!   sim_vz_pert   Scales [-1,1] random number in z direction: set to zero for uniform flow
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use Driver_interface, ONLY : Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none

#include "constants.h"
#include "Flash.h"
  
  
  integer :: i, j, status
 
  call Driver_getMype(MESH_COMM, sim_meshMe)

!! ambient values
  call RuntimeParameters_get("sim_p_amb",   sim_p_amb )
  call RuntimeParameters_get("sim_rho_amb", sim_rho_amb )
  call RuntimeParameters_get("sim_vx_amb",  sim_vx_amb )
  call RuntimeParameters_get("sim_vx_multiplier",  sim_vx_multiplier )

!! velocity perturbations
  call RuntimeParameters_get("sim_seed",   sim_seed )
  call RuntimeParameters_get("sim_vx_pert",   sim_vx_pert )
  call RuntimeParameters_get("sim_vy_pert",   sim_vy_pert )
  call RuntimeParameters_get("sim_vz_pert",   sim_vz_pert )
  call RuntimeParameters_get("sim_addPartCount",   sim_addpartCount )
  call RuntimeParameters_get("sim_addPartDisp",   sim_addpartDisp )
  call RuntimeParameters_get("gr_ptRemove",       sim_grPtRemove )
  call RuntimeParameters_get("xmin", sim_globalBndBox(LOW,IAXIS))
  call RuntimeParameters_get("xmax", sim_globalBndBox(HIGH,IAXIS))
  call RuntimeParameters_get("ymin", sim_globalBndBox(LOW,JAXIS))
  call RuntimeParameters_get("ymax", sim_globalBndBox(HIGH,JAXIS))
  call RuntimeParameters_get("zmin", sim_globalBndBox(LOW,KAXIS))
  call RuntimeParameters_get("zmax", sim_globalBndBox(HIGH,KAXIS))



  return
end subroutine Simulation_init
