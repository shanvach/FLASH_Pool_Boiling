!!****if* source/Simulation/SimulationMain/magnetoHD/LinearWave/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  Simulation_data()
!!
!!  Stores the local data for Simulation setup: Linearized MHD wave
!!  
!!  Reference: Crockett et al., JCP 203 (2005) 422-448!!
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"

  !! *** Runtime Parameters *** !!
  real, save    :: sim_gamma
  real, save    :: sim_smallP
  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  real, save    :: sim_nx, sim_ny, sim_delperturb
  real, save    :: sim_dens0,sim_B0,sim_pres0
  integer, save :: sim_choice
  logical, save :: sim_gCell, sim_killdivb, sim_steady
  character(len=MAX_STRING_LENGTH), save :: sim_choice_str

  integer, save :: sim_meshMe
end module Simulation_data
