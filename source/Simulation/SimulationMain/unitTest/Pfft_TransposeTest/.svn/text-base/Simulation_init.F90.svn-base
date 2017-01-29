!!****if* source/Simulation/SimulationMain/unitTest/Pfft_TransposeTest/Simulation_init
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
!!
!! ARGUMENTS
!!
!!   
!!
!! PARAMETERS
!!
!!  sim_outputGridData   Whether to create intermediate grid files.
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine Simulation_init()
  use Simulation_data, ONLY : sim_outputGridData, sim_meshMe
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_getMype
  implicit none
  call Driver_getMype(MESH_COMM,sim_meshMe)
  call RuntimeParameters_get('output_grid_data', sim_outputGridData)
end subroutine Simulation_init
