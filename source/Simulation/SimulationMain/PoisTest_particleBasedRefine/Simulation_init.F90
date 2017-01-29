!!****if* source/Simulation/SimulationMain/PoisTest_particleBasedRefine/Simulation_init
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
!!
!! DESCRIPTION
!!
!!  Initializes all the parameters needed for the Huang & Greengard (poistest)
!!  problem
!!
!! ARGUMENTS
!!
!!  
!!
!! PARAMETERS
!!
!!  sim_smlRho : 
!!  sim_gam : 
!!***

subroutine Simulation_init()
  
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp
  use Driver_interface, ONLY : Driver_getMype
  implicit none
#include "Flash.h"
#include "constants.h"
  

  call RuntimeParameters_get('sim_smlRho', sim_smlRho)
  
  call Logfile_stamp( "initializing for Huang & Greengard problem",  &
       "[Simulation_init]")

  sim_gam = 1.5
  call Driver_getMype(MESH_COMM, sim_meshMe)

end subroutine Simulation_init
