!!****if* source/Simulation/SimulationMain/RadShock/RadShock2d/Simulation_init
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
  use RuntimeParameters_interface, ONLY: RuntimeParameters_get
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use Simulation_data
  use Logfile_interface, ONLY : Logfile_stamp
  
  implicit none
#include "Flash.h"
#include "constants.h"

  call RuntimeParameters_get('sim_rhoXe' , sim_rhoXe )
  call RuntimeParameters_get('sim_teleXe', sim_teleXe)
  call RuntimeParameters_get('sim_tionXe', sim_tionXe)
  call RuntimeParameters_get('sim_tradXe', sim_tradXe) 
  call RuntimeParameters_get('sim_abarXe', sim_abarXe)
  call RuntimeParameters_get('sim_zbarXe', sim_zbarXe)

  call RuntimeParameters_get('sim_rhoCh' , sim_rhoCh )
  call RuntimeParameters_get('sim_teleCh', sim_teleCh)
  call RuntimeParameters_get('sim_tionCh', sim_tionCh)
  call RuntimeParameters_get('sim_tradCh', sim_tradCh)
  call RuntimeParameters_get('sim_abarCh', sim_abarCh)
  call RuntimeParameters_get('sim_zbarCh', sim_zbarCh)

  call RuntimeParameters_get('sim_rhoVa' , sim_rhoVa )
  call RuntimeParameters_get('sim_teleVa', sim_teleVa)
  call RuntimeParameters_get('sim_tionVa', sim_tionVa)
  call RuntimeParameters_get('sim_tradVa', sim_tradVa)
  call RuntimeParameters_get('sim_abarVa', sim_abarVa)
  call RuntimeParameters_get('sim_zbarVa', sim_zbarVa)

  call RuntimeParameters_get('sim_vely', sim_vely)
  call RuntimeParameters_get('sim_nbuffer', sim_nbuffer)
  call RuntimeParameters_get('rt_useMGD', sim_useMGD)
  call RuntimeParameters_get('sim_reflectDist', sim_reflectDist)

end subroutine Simulation_init
