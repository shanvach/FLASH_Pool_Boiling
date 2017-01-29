!!****if* source/Simulation/SimulationMain/Soundwave/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  Store the simulation data
!!
!! 
!!
!!
!! 
!!
!!
!!   
!!
!!***

module Simulation_data

  implicit none

  integer, save :: sim_meshMe

  !! *** Runtime Parameters *** !!

  real, save :: sim_gamma, sim_smallP, sim_smallX
  integer, save ::  sim_orientation
  real, save :: sim_wavelength, sim_rhoInit, sim_perturbAmp, sim_cs

  logical, save :: sim_gCell

end module Simulation_data


