!!****if* source/Simulation/SimulationMain/QuietPiston/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data for the QuietPiston problem
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!  sim_pAmbient    Initial ambient pressure
!!  sim_rhoAmbient  Initial ambient density
!!  sim_windVel     Inflow velocity (parallel to x-axis)
!!  gamma           the Gamma EOS thing
!!  smallp          minimum for pressure
!!  smallx          minimum for abundance
!!
!!
!!   
!!
!!***

module Simulation_data

  implicit none

  !! *** Runtime Parameters *** !!

  real, save :: sim_pAmbient, sim_rhoAmbient, sim_windVel
  real, save :: sim_gamma, sim_smallP, sim_smallX, sim_temp

integer, save :: sim_meshMe
end module Simulation_data


