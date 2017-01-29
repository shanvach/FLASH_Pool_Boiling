!!****if* source/Simulation/SimulationMain/magnetoHD/LinearWave/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_init()
!!
!! ARGUMENTS
!!
!!  none
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!
!!  Reference: Crockett et al., JCP 203 (2005) 422-448
!!
!!***

#define ALFVEN_WAVE  1
#define FAST_WAVE    2
#define SLOW_WAVE    3

subroutine Simulation_init()

  use Simulation_data
  use Driver_interface, ONLY : Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none

#include "constants.h"
#include "Flash.h"

  

  call Driver_getMype(MESH_COMM, sim_meshMe)

  call RuntimeParameters_get('gamma',   sim_gamma)
  call RuntimeParameters_get('xmin',    sim_xMin)
  call RuntimeParameters_get('ymin',    sim_yMin)
  call RuntimeParameters_get('zmin',    sim_zMin)
  call RuntimeParameters_get('xmax',    sim_xMax)
  call RuntimeParameters_get('ymax',    sim_yMax)
  call RuntimeParameters_get('zmax',    sim_zMax)

  call RuntimeParameters_get('nx',      sim_nx)
  call RuntimeParameters_get('ny',      sim_ny)
  call RuntimeParameters_get('steady',  sim_steady)
  call RuntimeParameters_get('dens0',   sim_dens0)
  call RuntimeParameters_get('B0',      sim_B0)
  call RuntimeParameters_get('pres0',   sim_pres0)
  call RuntimeParameters_get('choice',  sim_choice_str)
  if(trim(sim_choice_str) == "Alfven" .or. &
     trim(sim_choice_str) == "alfven") then
     sim_choice = ALFVEN_WAVE
  else if(trim(sim_choice_str) == "Fast" .or. &
          trim(sim_choice_str) == "fast" ) then
     sim_choice = FAST_WAVE
  else if (trim(sim_choice_str) == "Slow" .or. &
           trim(sim_choice_str) == "slow" ) then
     sim_choice = SLOW_WAVE
  else
     call Driver_abortFlash&
          ("[Simulation_init]: Unknown simulation case - 'Alfven','Fast', or 'Slow'")
  end if

  call RuntimeParameters_get('delperturb',sim_delperturb)

  call RuntimeParameters_get('killdivb', sim_killdivb)
  call RuntimeParameters_get('smallp',   sim_smallP)

  sim_gCell = .true.
  
end subroutine Simulation_init
