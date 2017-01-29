!!****if* source/Simulation/SimulationMain/unitTest/SolidBody/Simulation_init
!!
!! NAME
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!  Simulation_init( )
!!
!! DESCRIPTION
!!   Initialize all the runtime parameters needed for the Solid body unit test
!!
!! ARGUMENTS
!!
!!   
!!
!! PARAMETERS
!!
!!
!!***

subroutine Simulation_init()

  use Driver_interface, ONLY : Driver_getMyPE, Driver_getNumProcs
  use Simulation_data
  implicit none

#include "constants.h"
#include "Flash.h"
  
  call Driver_getMype(MESH_COMM, sim_meshMe)
  call Driver_getNumProcs(MESH_COMM, sim_numProcs)

end subroutine Simulation_init
