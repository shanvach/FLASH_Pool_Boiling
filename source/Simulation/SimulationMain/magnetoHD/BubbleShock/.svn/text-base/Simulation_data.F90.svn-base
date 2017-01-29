!!****if* source/Simulation/SimulationMain/magnetoHD/BubbleShock/Simulation_data
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
!!  Stores the local data for Simulation setup: OrszagTang
!!
!!***

module Simulation_data

  implicit none

#include "constants.h"

  !! *** Runtime Parameters *** !!
  real, save    :: sim_smallX, sim_smallRho
  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  real, save    :: sim_densT,sim_presT,sim_velxT,sim_velyT,sim_velzT,sim_magxT,sim_magyT,sim_magzT,&
                   sim_densB,sim_presB,sim_velxB,sim_velyB,sim_velzB,sim_magxB,sim_magyB,sim_magzB,&
                   sim_temp1,sim_temp2,sim_gammac1,sim_gammac2,sim_int1,sim_int2,&
                   sim_lposn,sim_bubbleRadius, sim_bubbleXCtr, sim_bubbleYCtr,   &
                   sim_bubbleZCtr,sim_bubbleDensity
  real, save    :: sim_gam1, sim_gam2
  logical, save :: sim_gCell, sim_killdivb

  integer, save :: sim_meshMe
end module Simulation_data
