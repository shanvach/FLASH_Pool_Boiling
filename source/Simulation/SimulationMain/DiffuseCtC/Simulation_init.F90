!!****if* source/Simulation/SimulationMain/DiffuseCtC/Simulation_init
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
#include "Flash.h"
  
  call RuntimeParameters_get('sim_rho' ,          sim_rho )

  call RuntimeParameters_get('sim_tele1',         sim_tele1)
  call RuntimeParameters_get('sim_tele2',         sim_tele2)
  call RuntimeParameters_get('sim_teleOffset',    sim_teleOffset)
  call RuntimeParameters_get('sim_teleSteepness', sim_teleSteepness)

  call RuntimeParameters_get('sim_tion1',         sim_tion1)
  call RuntimeParameters_get('sim_tion2',         sim_tion2)
  call RuntimeParameters_get('sim_tionOffset',    sim_tionOffset)
  call RuntimeParameters_get('sim_tionSteepness', sim_tionSteepness)

  call RuntimeParameters_get('sim_trad1',         sim_trad1)
  call RuntimeParameters_get('sim_trad2',         sim_trad2)
  call RuntimeParameters_get('sim_tradOffset',    sim_tradOffset)
  call RuntimeParameters_get('sim_tradSteepness', sim_tradSteepness)

  call RuntimeParameters_get('sim_initGeom', sim_initGeom)

  if(sim_initGeom == "polar" .and. NDIM == 1) then
     call Driver_abortFlash("Must use 2D or 3D for sim_initGeom == polar")
  end if

end subroutine Simulation_init
