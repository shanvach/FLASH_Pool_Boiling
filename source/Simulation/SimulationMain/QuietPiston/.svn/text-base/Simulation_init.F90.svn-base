!!****if* source/Simulation/SimulationMain/QuietPiston/Simulation_init
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
!!  Initializes all the parameters needed for the QuietPiston problem
!!
!! ARGUMENTS
!!
!!  
!!
!! PARAMETERS
!!
!!  sim_pAmbient    Initial ambient pressure
!!  sim_rhoAmbient  Initial ambient density
!!  sim_windVel     Inflow velocity (parallel to x-axis)
!!  gamma           the Gamma EOS thing
!!  smallp          minimum for pressure
!!  smallx          minimum for abundance
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp
  implicit none

#include "constants.h"
#include "Flash.h"

  

  call RuntimeParameters_get('sim_pAmbient', sim_pAmbient)
  call RuntimeParameters_get('sim_rhoAmbient', sim_rhoAmbient)
  call RuntimeParameters_get('sim_windVel', sim_windVel)
  call RuntimeParameters_get('sim_temp', sim_temp)
  call RuntimeParameters_get('gamma', sim_gamma)
  call RuntimeParameters_get('smallx', sim_smallX)
  call RuntimeParameters_get('smallp', sim_smallP)


  call Logfile_stamp("initializing for QuietPiston", 'run_init')
!!  write (*,*) "flash:  initializing for QuietPiston"


end subroutine Simulation_init
