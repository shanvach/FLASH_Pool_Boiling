!!****if* source/Simulation/SimulationMain/magnetoHD/RTmhd/Simulation_init
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
!!  Initializes initial conditions for MHD KH instability problem.
!!
!!***

subroutine Simulation_init()

  use Simulation_data, ONLY : sim_gamma, sim_xMin, sim_yMin, sim_zMin, &
                              sim_xMax, sim_yMax, sim_zMax, sim_gCell, &
                              sim_killdivb, sim_Bx0, sim_By0, sim_Bz0, &
                              sim_epsilon,sim_gconst, sim_rhoHeavy, sim_rhoLight,&
                              sim_number, sim_smallP, sim_smallX,sim_meshMe

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_getMype

  implicit none

#include "constants.h"
#include "Flash.h"

  

  call RuntimeParameters_get('gamma',   sim_gamma)
  call RuntimeParameters_get('xmin',    sim_xMin)
  call RuntimeParameters_get('ymin',    sim_yMin)
  call RuntimeParameters_get('zmin',    sim_zMin)
  call RuntimeParameters_get('xmax',    sim_xMax)
  call RuntimeParameters_get('ymax',    sim_yMax)
  call RuntimeParameters_get('zmax',    sim_zMax)
  call RuntimeParameters_get('gconst', sim_gconst)
  call RuntimeParameters_get('epsilon',sim_epsilon)
  call RuntimeParameters_get('rho_heavy',sim_rhoHeavy)
  call RuntimeParameters_get('rho_light',sim_rhoLight)
  call RuntimeParameters_get('simulation',sim_number)

#if defined(MAGX_VAR) && defined(MAGY_VAR) && defined(MAGZ_VAR) && defined(DIVB_VAR)
  call RuntimeParameters_get('killdivb',sim_killdivb)
  call RuntimeParameters_get('smallp', sim_smallP)
  call RuntimeParameters_get('smallx', sim_smallX)
  call RuntimeParameters_get('Bx0',    sim_Bx0)
  call RuntimeParameters_get('By0',    sim_By0)
  call RuntimeParameters_get('Bz0',    sim_Bz0)
#endif
  call Driver_getMype(MESH_COMM, sim_meshMe)
  sim_gCell = .true.

end subroutine Simulation_init
