!!****if* source/Simulation/SimulationMain/BuriedCone/Simulation_init
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
  use Logfile_interface, ONLY : Logfile_stamp
  
  implicit none

#include "constants.h"
#include "Flash.h"

  real :: xmin, xmax, ymin, ymax
  integer :: lrefine_max, nblockx, nblocky
  character(len=MAX_STRING_LENGTH) :: str



  call RuntimeParameters_get('sim_targetRadius', sim_targetRadius)
  call RuntimeParameters_get('sim_targetHeight', sim_targetHeight)
  call RuntimeParameters_get('sim_vacuumHeight', sim_vacuumHeight)
  call RuntimeParameters_get('sim_coneTipRadius', sim_coneTipRadius)
  call RuntimeParameters_get('sim_coneDepth', sim_coneDepth)
  call RuntimeParameters_get('sim_coneStart', sim_coneStart)
  
  call RuntimeParameters_get('sim_rhoTarg', sim_rhoTarg)
  call RuntimeParameters_get('sim_teleTarg', sim_teleTarg)
  call RuntimeParameters_get('sim_tionTarg', sim_tionTarg)
  call RuntimeParameters_get('sim_tradTarg', sim_tradTarg)
  call RuntimeParameters_get('ms_targA', sim_abarTarg)
  call RuntimeParameters_get('ms_targZ', sim_zbarTarg)
  
  call RuntimeParameters_get('sim_rhoCham', sim_rhoCham)
  call RuntimeParameters_get('sim_teleCham', sim_teleCham)
  call RuntimeParameters_get('sim_tionCham', sim_tionCham)
  call RuntimeParameters_get('sim_tradCham', sim_tradCham)
  call RuntimeParameters_get('ms_chamA', sim_abarCham)
  call RuntimeParameters_get('ms_chamZ', sim_zbarCham)

  call RuntimeParameters_get('smallX', sim_smallX)

end subroutine Simulation_init
