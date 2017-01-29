!!****if* source/Simulation/SimulationMain/magnetoHD/AdvectBall/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  Simulation_data()
!!
!!  Stores the local data for Simulation setup: FieldLoop
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
  real, save    :: sim_xCtr, sim_yCtr, sim_zCtr
  real, save    :: sim_rx, sim_ry
  real, save    :: sim_Ux_initial,sim_Uy_initial,sim_Uz_initial, sim_ballRadius

  logical, save :: sim_gCell, sim_killdivb

  integer, save :: sim_meshMe
end module Simulation_data
