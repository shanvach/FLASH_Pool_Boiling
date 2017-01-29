!!****if* source/Simulation/SimulationMain/magnetoHD/KHmhd/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_init()
!!
!! ARGUMENTS
!!
!!  none
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Initializes initial conditions for MHD KH instability problem.
!!
!!***

subroutine Simulation_init()

  use Simulation_data, ONLY : sim_gamma, sim_xMin, sim_yMin, sim_zMin, &
                              sim_xMax, sim_yMax, sim_zMax, sim_gCell, &
                              sim_killdivb, sim_Vx0, sim_Bx0, sim_By0, &
                              sim_Bz0, sim_epsilon, sim_smallP, sim_smallX,&
                              sim_meshMe

  use Driver_interface, ONLY : Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "constants.h"
#include "Flash.h"

  
  call Driver_getMype(MESH_COMM, sim_meshMe)

  call RuntimeParameters_get('gamma',   sim_gamma)
  call RuntimeParameters_get('xmin',    sim_xMin)
  call RuntimeParameters_get('ymin',    sim_yMin)
  call RuntimeParameters_get('zmin',    sim_zMin)
  call RuntimeParameters_get('xmax',    sim_xMax)
  call RuntimeParameters_get('ymax',    sim_yMax)
  call RuntimeParameters_get('zmax',    sim_zMax)
  call RuntimeParameters_get('killdivb',sim_killdivb)
  call RuntimeParameters_get('smallp',  sim_smallP)
  call RuntimeParameters_get('smallx',  sim_smallX)

  call RuntimeParameters_get('Vx0',    sim_Vx0)
  call RuntimeParameters_get('Bx0',    sim_Bx0)
  call RuntimeParameters_get('By0',    sim_By0)
  call RuntimeParameters_get('Bz0',    sim_Bz0)
  call RuntimeParameters_get('epsilon',sim_epsilon)

  sim_gCell = .true.

end subroutine Simulation_init
