!!****if* source/Simulation/SimulationMain/unitTest/Gravity/Poisson2/Simulation_init
!!
!! NAME
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!  Simulation_init( )
!!
!! DESCRIPTION   
!!   Initialize all the runtime parameters needed for the Particle unitTest
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
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  implicit none

#include "constants.h"
#include "Flash.h"
  
  
 

  !-------------------------------------------------------------------------
  !               Write a message to stdout describing the problem setup.
  !-------------------------------------------------------------------------
  


  
  call RuntimeParameters_get('discRadius',sim_discRadius)
  call RuntimeParameters_get('density',sim_density)
  call PhysicalConstants_get("Newton",sim_newton)
  return
end subroutine Simulation_init
