!!****if* source/Simulation/SimulationMain/SamraiTest/Simulation_init
!!
!! NAME
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!  call Simulation_init()
!!
!!
!! DESCRIPTION
!!  Initializes all the parameters needed for a particular simulation
!!
!! ARGUMENTS
!!
!!  none 
!!
!!***

subroutine Simulation_init()
  
  use Simulation_data 
  use Driver_interface, ONLY : Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Grid_interface, ONLY:  Grid_getNumProcs

  implicit none

#include "Flash.h"
#include "constants.h"

  ! Everybody should know these
  call Driver_getMype(MESH_COMM, sim_meshMe)
  call Grid_getNumProcs(sim_numProcs)

#ifdef FIXEDBLOCKSIZE
  sim_gridIMax = GRID_IHI_GC
  sim_gridJMax = GRID_JHI_GC
  sim_gridKMax = GRID_KHI_GC
  call RuntimeParameters_get('iGridSize',sim_iGridSize)
  call RuntimeParameters_get('jGridSize',sim_jGridSize)
  call RuntimeParameters_get('kGridSize',sim_kGridSize)

#else
  call RuntimeParameters_get('gridImax',sim_gridIMax)
  call RuntimeParameters_get('gridJmax',sim_gridJMax)
  call RuntimeParameters_get('gridKmax',sim_gridKMax)
  call RuntimeParameters_get('iGridSize',sim_iGridSize)
  call RuntimeParameters_get('jGridSize',sim_jGridSize)
  call RuntimeParameters_get('kGridSize',sim_kGridSize)
#endif
  sim_gCell = .false.

end subroutine Simulation_init







