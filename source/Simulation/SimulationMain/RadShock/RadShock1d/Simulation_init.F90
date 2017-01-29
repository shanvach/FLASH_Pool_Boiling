!!****if* source/Simulation/SimulationMain/RadShock/RadShock1d/Simulation_init
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
  use RadTrans_interface, ONLY : RadTrans_mgdSetBc
  
  implicit none

#include "constants.h"

  call RuntimeParameters_get('sim_rho' , sim_rho )
  call RuntimeParameters_get('sim_tele', sim_tele)
  call RuntimeParameters_get('sim_tion', sim_tion)
  call RuntimeParameters_get('sim_trad', sim_trad)
  call RuntimeParameters_get('sim_velx', sim_velx)
  call RuntimeParameters_get('smallX',   sim_smallX)

  call RuntimeParameters_get('sim_slabThickness', sim_slabThickness)
  call RuntimeParameters_get('sim_rhoBe' , sim_rhoBe )
  call RuntimeParameters_get('sim_teleBe', sim_teleBe)
  call RuntimeParameters_get('sim_tionBe', sim_tionBe)
  call RuntimeParameters_get('sim_tradBe', sim_tradBe)

  call RuntimeParameters_get('sim_vacThickness', sim_vacThickness)
  call RuntimeParameters_get('sim_rhoVa' , sim_rhoVa )
  call RuntimeParameters_get('sim_teleVa', sim_teleVa)
  call RuntimeParameters_get('sim_tionVa', sim_tionVa)
  call RuntimeParameters_get('sim_tradVa', sim_tradVa)

  call RuntimeParameters_get('sim_specialGroup', sim_specialGroup)
  call RuntimeParameters_get('sim_specialUrad', sim_specialUrad)

  if(sim_specialGroup > 0) then
     call RadTrans_mgdSetBc(ig=sim_specialGroup, f=1, &
          bcType=DIRICHLET, bcValue=sim_specialUrad)
  end if

end subroutine Simulation_init
