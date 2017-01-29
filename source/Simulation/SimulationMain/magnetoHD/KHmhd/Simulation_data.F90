!!****if* source/Simulation/SimulationMain/magnetoHD/KHmhd/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data
!!
!! DESCRIPTION
!!
!!  Stores the local data for Simulation setup: MHD KH
!!
!!***

module Simulation_data

  implicit none

#include "constants.h"

  !! *** Runtime Parameters *** !!
  real, save    :: sim_gamma,sim_smallX, sim_smallRho, sim_smallP
  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  real, save    :: sim_Vx0, sim_Bx0, sim_By0, sim_Bz0, sim_epsilon
  logical, save :: sim_gCell, sim_killdivb

  integer, save :: sim_meshMe
end module Simulation_data
