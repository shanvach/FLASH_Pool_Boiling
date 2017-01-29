!!****if* source/Simulation/SimulationMain/magnetoHD/BierSod/Simulation_init
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
!!  Initializes all the parameters needed for a particular simulation
!!
!!
!! ARGUMENTS
!!
!!  
!!
!! PARAMETERS
!!
!!***

subroutine Simulation_init()
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none

#include "constants.h"
#include "Flash.h"
  
  call RuntimeParameters_get('sim_densTarg', sim_densTarg)
  call RuntimeParameters_get('sim_teleTarg', sim_teleTarg)
  call RuntimeParameters_get('sim_tionTarg', sim_tionTarg)
  
  call RuntimeParameters_get('sim_densCham', sim_densCham)
  call RuntimeParameters_get('sim_teleCham', sim_teleCham)
  call RuntimeParameters_get('sim_tionCham', sim_tionCham)

  call RuntimeParameters_get('smallX', sim_smallX)
  call RuntimeParameters_get('sim_ndiv', sim_ndiv)
  call RuntimeParameters_get('sim_radius', sim_radius)

end subroutine Simulation_init
