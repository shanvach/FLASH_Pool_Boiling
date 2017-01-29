!!****if* source/Simulation/SimulationMain/magnetoHD/CloudFieldLine/Simulation_data
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
!!  Stores the local data for Simulation setup: BrioWu
!!
!!***

module Simulation_data

  implicit none

#include "constants.h"

  !! *** Runtime Parameters *** !!    
  real, save   :: sim_posn, sim_possh,sim_gamma
  real, save   :: sim_uRight, sim_uLeft, sim_vRight, sim_vLeft, sim_wLeft, sim_wRight
  real, save   :: sim_dRight, sim_dLeft, sim_pRight, sim_pLeft
  real, save   :: sim_xmin, sim_xmax, sim_ymin, sim_ymax
  real, save   :: sim_cloudRadius, sim_cloudXCtr, sim_cloudYCtr,   &
                  sim_cloudZCtr,sim_cloudDensity

  !! Simulation variables
  real, save    :: sim_xangle, sim_yangle, sim_xcos, sim_ycos, sim_zcos
  real, save    :: sim_smallx, sim_smallP
  real, save    :: sim_rx, sim_ry
  logical, save :: sim_gCell, sim_killdivb

  integer, save :: sim_meshMe
end module Simulation_data
