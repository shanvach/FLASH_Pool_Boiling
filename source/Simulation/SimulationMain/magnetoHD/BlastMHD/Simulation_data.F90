!!****if* source/Simulation/SimulationMain/magnetoHD/BlastMHD/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  Simulation_data()
!!
!!  Stores the local data for Simulation setup: MHD Blast wave problem
!!  
!!  Reference:  Gardiner & Stone JCP 205(2005),509-539
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"

  !! *** Variables pertaining to this Simulation *** !!
  real, save    :: sim_gamma
  real, save    :: sim_radius, sim_beta_inner, sim_beta_outer
  real, save    :: sim_velx, sim_vely, sim_velz
  real, save    :: sim_magx, sim_magy, sim_magz
  real, save    :: sim_smallP
  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  logical, save :: sim_gCell, sim_killdivb

  integer, save :: sim_meshMe
end module Simulation_data
