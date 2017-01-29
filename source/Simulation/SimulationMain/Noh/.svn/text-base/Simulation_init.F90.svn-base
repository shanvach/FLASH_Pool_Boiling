!!****if* source/Simulation/SimulationMain/Noh/Simulation_init
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
!!  Initializes all the parameters needed for the Sod shock tube
!!  problem
!!
!! ARGUMENTS
!!
!!  
!!
!! PARAMETERS
!!
!!
!!***

subroutine Simulation_init()
  
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp
  implicit none
!!#include "RuntimeParameters.h"
#include "Flash.h"
!!#include "Logfile.h"

  

  call RuntimeParameters_get('sim_rhoInit', sim_rhoInit)  
  call RuntimeParameters_get('sim_pInit', sim_pInit)
  call RuntimeParameters_get('sim_uInit', sim_uInit)
  call RuntimeParameters_get('sim_gamma', sim_gamma)

  call RuntimeParameters_get('xmin', xmin)
  call RuntimeParameters_get('xmax', xmax)
  call RuntimeParameters_get('ymin', ymin)
  call RuntimeParameters_get('ymax', ymax)
  call RuntimeParameters_get('zmin', zmin)
  call RuntimeParameters_get('zmax', zmax)

  call RuntimeParameters_get('smallp', smallp)

  call Logfile_stamp( "initializing Noh problem",  &
       "[Simulation_init]")
     
end subroutine Simulation_init







