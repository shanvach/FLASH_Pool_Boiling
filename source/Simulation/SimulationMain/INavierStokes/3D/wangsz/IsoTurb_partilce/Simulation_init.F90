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
                              sim_xMax, sim_yMax, sim_zMax, sim_gCell, &
                              sim_fbao, sim_turbkin_expect, sim_scaleForce

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use Grid_data, only : gr_meshMe

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

  call RuntimeParameters_get('sim_fbao', sim_fbao)
  call RuntimeParameters_get('sim_turbkin_expect', sim_turbkin_expect)
  call RuntimeParameters_get('sim_scaleForce', sim_scaleForce)
  if(gr_meshMe == MASTER_PE) write(*,*) 'Forced turbulence, FBAO=', sim_fbao
  if(gr_meshMe == MASTER_PE) write(*,*) 'Expected tke =', sim_turbkin_expect
  if(gr_meshMe == MASTER_PE) write(*,*) 'Initial force scalor =', sim_scaleForce

end subroutine Simulation_init
