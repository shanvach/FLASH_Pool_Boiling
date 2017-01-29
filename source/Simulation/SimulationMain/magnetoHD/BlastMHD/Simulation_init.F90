!!****if* source/Simulation/SimulationMain/magnetoHD/BlastMHD/Simulation_init
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
!!   none
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!
!!  Reference:  Gardiner & Stone JCP 205(2005),509-539
!!
!!***

subroutine Simulation_init()

  use Simulation_data
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

  call RuntimeParameters_get('killdivb',   sim_killdivb)
  call RuntimeParameters_get('smallp',     sim_smallP)
  call RuntimeParameters_get('radius',     sim_radius)
  call RuntimeParameters_get('beta_inner', sim_beta_inner)
  call RuntimeParameters_get('beta_outer', sim_beta_outer)
  call RuntimeParameters_get('velx_init',  sim_velx)
  call RuntimeParameters_get('vely_init',  sim_vely)
  call RuntimeParameters_get('velz_init',  sim_velz)
  call RuntimeParameters_get('magx_init',  sim_magx)
  call RuntimeParameters_get('magy_init',  sim_magy)
  call RuntimeParameters_get('magz_init',  sim_magz)

  sim_gCell = .true.
  
end subroutine Simulation_init
