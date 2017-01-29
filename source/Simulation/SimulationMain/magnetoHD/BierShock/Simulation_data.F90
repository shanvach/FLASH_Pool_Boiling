!!****if* source/Simulation/SimulationMain/magnetoHD/BierShock/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!  Use Simulation_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data
!!
!! 
!!***

module Simulation_data

  implicit none

  !! *** Simulation Runtime Parameters *** !!
  real, save :: sim_skewFactor
  real, save :: sim_rho
  real, save :: sim_eint
  real, save :: sim_velx
  real, save :: sim_erad

  !! *** Other Runtime Parameters *** !!
  real, save :: sim_singleSpeciesZ
  real, save :: sim_ymin
  real, save :: sim_ymax

end module Simulation_data


