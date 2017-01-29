!!****if* source/Simulation/SimulationMain/magnetoHD/CloudShock/Simulation_data
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
!!  Stores the local data for Simulation setup: CloudShock
!!
!!***

module Simulation_data

  implicit none

#include "constants.h"

  !! *** Runtime Parameters *** !!
  real, save    :: sim_gamma,sim_smallX, sim_smallRho, sim_smallP
  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  real, save    :: sim_dR,sim_pR,sim_uR,sim_vR,sim_wR,sim_bxR,sim_byR,sim_bzR,&
                   sim_dL,sim_pL,sim_uL,sim_vL,sim_wL,sim_bxL,sim_byL,sim_bzL,&
                   sim_lposn,sim_cloudRadius, sim_cloudXCtr, sim_cloudYCtr,   &
                   sim_cloudZCtr,sim_cloudDensity
  logical, save :: sim_gCell, sim_killdivb

  integer, save :: sim_meshMe, sim_Type
end module Simulation_data
