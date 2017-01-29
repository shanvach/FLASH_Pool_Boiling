!!****if* source/Simulation/SimulationMain/unitTest/Particles/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!  module Simulation_data()
!!
!! DESCRIPTION
!!
!!  Store the simulation data for unitTesting of Particles
!!   
!! ARGUMENTS
!!
!! PARAMETERS   
!!      Described in the Config file
!!
!! NOTES
!!
!!  No arguments.  All data passed by "use Simulation_data"
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"
#include "Flash.h"
  
  real, save    :: sim_rho_amb, sim_p_amb, sim_vx_amb, sim_vx_multiplier
  real, save    :: sim_seed, sim_vx_pert, sim_vy_pert, sim_vz_pert, sim_addPartDisp
  real, dimension(LOW:HIGH,MDIM), save    :: sim_globalBndBox

  integer, save :: sim_meshMe, sim_addPartCount

  logical, save :: sim_grPtRemove !! for gr_ptRemove RP
end module Simulation_data
