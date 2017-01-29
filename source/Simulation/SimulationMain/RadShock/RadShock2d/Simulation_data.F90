!!****if* source/Simulation/SimulationMain/RadShock/RadShock2d/Simulation_data
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

  !! *** Runtime Parameters *** !!
  integer, save :: sim_nbuffer
  real, save :: sim_reflectDist
  real, save :: sim_vely
  logical, save :: sim_useMGD
  real, save :: sim_radconst

  real, save :: sim_rhoXe
  real, save :: sim_tionXe
  real, save :: sim_teleXe
  real, save :: sim_tradXe
  real, save :: sim_zbarXe
  real, save :: sim_abarXe

  real, save :: sim_rhoCh
  real, save :: sim_tionCh
  real, save :: sim_teleCh
  real, save :: sim_tradCh
  real, save :: sim_zbarCh
  real, save :: sim_abarCh

  real, save :: sim_rhoVa
  real, save :: sim_tionVa
  real, save :: sim_teleVa
  real, save :: sim_tradVa
  real, save :: sim_zbarVa
  real, save :: sim_abarVa

end module Simulation_data


