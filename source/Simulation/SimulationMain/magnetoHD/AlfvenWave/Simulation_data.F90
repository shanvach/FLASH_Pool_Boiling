!!****if* source/Simulation/SimulationMain/magnetoHD/AlfvenWave/Simulation_data
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
!!  Stores the local data for Simulation setup: Circularly polarized Alfven wave
!!  
!!  Reference: Gardiner & Stone JCP 205(2005),509-539
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"

  !! *** Runtime Parameters *** !!
  real, save    :: sim_gamma
  real, save    :: sim_smallRho, sim_smallP
  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  real, save    :: sim_rx, sim_ry
  real, save    :: sim_U0,sim_B0,sim_P0

  logical, save :: sim_gCell, sim_killdivb, sim_steady

  integer, save :: sim_meshMe
end module Simulation_data
