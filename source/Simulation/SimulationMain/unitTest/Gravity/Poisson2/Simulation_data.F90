!!****if* source/Simulation/SimulationMain/unitTest/Gravity/Poisson2/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!  module Simulation_data()
!!
!! DESCRIPTION
!!
!!  Store the simulation data for unitTesting of Poisson solver
!!   
!! ARGUMENTS
!!
!! PARAMETERS   
!!      Described in the Config file
!!
!! NOTES
!!
!!  No arguments.  All data passed by "use Simulation_data"
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"
#include "Flash.h"
  
  real, save :: sim_discRadius,sim_density,sim_newton

integer, save :: sim_meshMe
end module Simulation_data
