!!****if* source/Simulation/SimulationMain/magnetoHD/SupersonicJets/Simulation_init
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
!!  Initializes initial conditions for Supersonic Jet problem.
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
  call RuntimeParameters_get('killdivb',sim_killdivb)
  call RuntimeParameters_get('smallp',  sim_smallP)
  call RuntimeParameters_get('smallx',  sim_smallX)

  call RuntimeParameters_get('d_ambient',sim_dAmbient)
  call RuntimeParameters_get('p_ambient',sim_pAmbient)
  call RuntimeParameters_get('u_ambient',sim_uAmbient)
  call RuntimeParameters_get('v_ambient',sim_vAmbient)
  call RuntimeParameters_get('w_ambient',sim_wAmbient)
  call RuntimeParameters_get('bx_ambient',sim_bxAmbient)
  call RuntimeParameters_get('by_ambient',sim_byAmbient)
  call RuntimeParameters_get('bz_ambient',sim_bzAmbient)

  call RuntimeParameters_get('jetRadius',sim_jetRadius)
  call RuntimeParameters_get('jetXCtr',sim_jetXCtr)
  call RuntimeParameters_get('jetYCtr',sim_jetYCtr)
  call RuntimeParameters_get('jetZCtr',sim_jetZCtr)
  call RuntimeParameters_get('jetDensity',sim_jetDensity)
  call RuntimeParameters_get('jetVelocity',sim_jetVelocity)

  sim_gCell = .true.

end subroutine Simulation_init
