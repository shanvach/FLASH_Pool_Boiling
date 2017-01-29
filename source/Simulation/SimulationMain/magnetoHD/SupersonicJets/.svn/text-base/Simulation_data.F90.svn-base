!!****if* source/Simulation/SimulationMain/magnetoHD/SupersonicJets/Simulation_data
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
!!  Stores the local data for Simulation setup: OrszagTang
!!
!!***

module Simulation_data

  implicit none

#include "constants.h"

  integer, save :: sim_meshMe

  !! *** Runtime Parameters *** !!
  real, save    :: sim_gamma,sim_smallX, sim_smallRho, sim_smallP
  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  real, save    :: sim_dAmbient,sim_pAmbient,sim_uAmbient,sim_vAmbient,sim_wAmbient,&
                   sim_bxAmbient,sim_byAmbient,sim_bzAmbient,&
                   sim_lposn,sim_jetRadius, sim_jetXCtr, sim_jetYCtr,   &
                   sim_jetZCtr,sim_jetDensity, sim_jetVelocity
  logical, save :: sim_gCell, sim_killdivb
end module Simulation_data
