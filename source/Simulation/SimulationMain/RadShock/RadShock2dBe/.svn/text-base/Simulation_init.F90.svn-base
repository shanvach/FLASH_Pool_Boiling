!!****if* source/Simulation/SimulationMain/RadShock/RadShock2dBe/Simulation_init
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
  use Simulation_data
  use Logfile_interface, ONLY : Logfile_stamp

#ifdef MGD_NGROUPS
  use RadTrans_interface, ONLY: RadTrans_mgdSetBound
#endif
  
  implicit none
#include "Flash.h"
#include "constants.h"

  
  integer :: g
  real, parameter :: K = 1.60217653e-12 ! Boltzmann Constant [erg/eV]

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

  call RuntimeParameters_get('sim_rhoBe' , sim_rhoBe )
  call RuntimeParameters_get('sim_teleBe', sim_teleBe)
  call RuntimeParameters_get('sim_tionBe', sim_tionBe)
  call RuntimeParameters_get('sim_tradBe', sim_tradBe)
  call RuntimeParameters_get('sim_abarBe', sim_abarBe)
  call RuntimeParameters_get('sim_zbarBe', sim_zbarBe)

  call RuntimeParameters_get('sim_gradSize', sim_gradSize)

  call RuntimeParameters_get('sim_tubeRadius', sim_tubeRadius)
  call RuntimeParameters_get('sim_tubeThickness', sim_tubeThickness)
  call RuntimeParameters_get('sim_slabThickness', sim_slabThickness)
  call RuntimeParameters_get('sim_slabRadius', sim_slabRadius)
  call RuntimeParameters_get('sim_vacThickness', sim_vacThickness)

end subroutine Simulation_init
