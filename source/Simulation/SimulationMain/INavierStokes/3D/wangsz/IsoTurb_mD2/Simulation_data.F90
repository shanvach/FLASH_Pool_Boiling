!!****if* source/Simulation/SimulationMain/INavierStokes/3D/ChannelLam/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  Simulation_data()
!!
!! DESCRIPTION
!!
!!  Stores the local data for Simulation setup: INS-iso-turb
!!
!!***

module Simulation_data

  implicit none

#include "constants.h"

  !! *** Runtime Parameters *** !!
  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  logical, save :: sim_gCell
  
  integer, save :: sim_meshMe
  real, save :: sim_turbkin_expect, sim_fbao, sim_scaleForce

end module Simulation_data
