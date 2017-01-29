!!****if* source/Simulation/SimulationMain/INavierStokes/3D/ChannelLam/Simulation_init
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
!! ARGUMENTS
!!
!!   
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Initializes initial conditions for INS-isotropic turbulence problem.
!!
!!***

subroutine Simulation_init()

  use Simulation_data, ONLY : sim_xMin, sim_yMin, sim_zMin, &
                              sim_xMax, sim_yMax, sim_zMax, sim_gCell

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "constants.h"
#include "Flash.h"

  

  call RuntimeParameters_get('xmin',    sim_xMin)
  call RuntimeParameters_get('ymin',    sim_yMin)
  call RuntimeParameters_get('zmin',    sim_zMin)
  call RuntimeParameters_get('xmax',    sim_xMax)
  call RuntimeParameters_get('ymax',    sim_yMax)
  call RuntimeParameters_get('zmax',    sim_zMax)

  sim_gCell = .true.

end subroutine Simulation_init
