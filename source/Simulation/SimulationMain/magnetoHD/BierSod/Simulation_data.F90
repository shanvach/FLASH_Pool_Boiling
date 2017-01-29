!!****if* source/Simulation/SimulationMain/magnetoHD/BierSod/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!  Use Simulation_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data
!!
!! 
!!***
module Simulation_data

  implicit none

#include "constants.h"

  !! *** Runtime Parameters *** !!    
  real, save :: sim_densTarg
  real, save :: sim_teleTarg
  real, save :: sim_tionTarg
  
  real, save :: sim_densCham
  real, save :: sim_teleCham
  real, save :: sim_tionCham

  real, save :: sim_smallX
  real, save :: sim_radius
  
  integer, save :: sim_ndiv

end module Simulation_data


