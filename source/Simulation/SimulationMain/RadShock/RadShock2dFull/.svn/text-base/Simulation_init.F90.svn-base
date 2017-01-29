!!****if* source/Simulation/SimulationMain/RadShock/RadShock2dFull/Simulation_init
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
  
  implicit none
#include "Flash.h"
#include "constants.h"

  
  integer :: g
  real, parameter :: K = 1.60217653e-12 ! Boltzmann Constant [erg/eV]

  call RuntimeParameters_get('smallx' , sim_smallX)

  call RuntimeParameters_get('sim_rhoXe' , sim_rhoXe )
  call RuntimeParameters_get('sim_teleXe', sim_teleXe)
  call RuntimeParameters_get('sim_tionXe', sim_tionXe)
  call RuntimeParameters_get('sim_tradXe', sim_tradXe) 
  call RuntimeParameters_get('ms_xeA', sim_abarXe)
  call RuntimeParameters_get('ms_xeZ', sim_zbarXe)

  call RuntimeParameters_get('sim_rhoCh' , sim_rhoCh )
  call RuntimeParameters_get('sim_teleCh', sim_teleCh)
  call RuntimeParameters_get('sim_tionCh', sim_tionCh)
  call RuntimeParameters_get('sim_tradCh', sim_tradCh)
  call RuntimeParameters_get('ms_poliA', sim_abarCh)
  call RuntimeParameters_get('ms_poliZ', sim_zbarCh)

  call RuntimeParameters_get('sim_rhoVa' , sim_rhoVa )
  call RuntimeParameters_get('sim_teleVa', sim_teleVa)
  call RuntimeParameters_get('sim_tionVa', sim_tionVa)
  call RuntimeParameters_get('sim_tradVa', sim_tradVa)
  call RuntimeParameters_get('ms_vacuA', sim_abarVa)
  call RuntimeParameters_get('ms_vacuZ', sim_zbarVa)

  call RuntimeParameters_get('sim_rhoBe' , sim_rhoBe )
  call RuntimeParameters_get('sim_teleBe', sim_teleBe)
  call RuntimeParameters_get('sim_tionBe', sim_tionBe)
  call RuntimeParameters_get('sim_tradBe', sim_tradBe)
  call RuntimeParameters_get('ms_beA', sim_abarBe)
  call RuntimeParameters_get('ms_beZ', sim_zbarBe)

  call RuntimeParameters_get('sim_rhoAu' , sim_rhoAu )
  call RuntimeParameters_get('sim_teleAu', sim_teleAu)
  call RuntimeParameters_get('sim_tionAu', sim_tionAu)
  call RuntimeParameters_get('sim_tradAu', sim_tradAu)
  call RuntimeParameters_get('ms_goldA', sim_abarAu)
  call RuntimeParameters_get('ms_goldZ', sim_zbarAu)

  call RuntimeParameters_get('sim_rhoAc' , sim_rhoAc )
  call RuntimeParameters_get('sim_teleAc', sim_teleAc)
  call RuntimeParameters_get('sim_tionAc', sim_tionAc)
  call RuntimeParameters_get('sim_tradAc', sim_tradAc)
  call RuntimeParameters_get('ms_acryA', sim_abarAc)
  call RuntimeParameters_get('ms_acryZ', sim_zbarAc)

  call RuntimeParameters_get('sim_tubeRadius', sim_tubeRadius)
  call RuntimeParameters_get('sim_tubeThickness', sim_tubeThickness)
  call RuntimeParameters_get('sim_slabThickness', sim_slabThickness)
  call RuntimeParameters_get('sim_slabRadius', sim_slabRadius)
  call RuntimeParameters_get('sim_vacThickness', sim_vacThickness)
  call RuntimeParameters_get('sim_goldThickness', sim_goldThickness)
  call RuntimeParameters_get('sim_goldRadius', sim_goldRadius)
  call RuntimeParameters_get('sim_acrylicThickness', sim_acrylicThickness)
  call RuntimeParameters_get('sim_acrylicRadius', sim_acrylicRadius)
  call RuntimeParameters_get('sim_windowThickness', sim_windowThickness)

  call RuntimeParameters_get('sim_lrefmaxBase', sim_lrefmaxBase)  

  call RuntimeParameters_get('sim_lrefmaxBe', sim_lrefmaxBe)
  call RuntimeParameters_get('sim_belrXMin' , sim_belrXMin )
  call RuntimeParameters_get('sim_belrXMax' , sim_belrXMax )
  call RuntimeParameters_get('sim_belrYMin' , sim_belrYMin )
  call RuntimeParameters_get('sim_belrYMax' , sim_belrYMax )
  
  call RuntimeParameters_get('sim_lrefmaxPoly', sim_lrefmaxPoly)
  call RuntimeParameters_get('sim_polylrXMin' , sim_polylrXMin )
  call RuntimeParameters_get('sim_polylrXMax' , sim_polylrXMax )
  call RuntimeParameters_get('sim_polylrYMin' , sim_polylrYMin )
  call RuntimeParameters_get('sim_polylrYMax' , sim_polylrYMax )

end subroutine Simulation_init
