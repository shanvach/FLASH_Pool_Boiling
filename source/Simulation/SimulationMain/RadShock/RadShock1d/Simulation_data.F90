!!****if* source/Simulation/SimulationMain/RadShock/RadShock1d/Simulation_data
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
  real, save :: sim_rho
  real, save :: sim_tion
  real, save :: sim_tele
  real, save :: sim_trad
  real, save :: sim_velx

  real, save :: sim_smallX

  integer, save :: sim_specialGroup
  real, save :: sim_specialUrad

  real, save :: sim_slabThickness
  real, save :: sim_rhoBe
  real, save :: sim_teleBe
  real, save :: sim_tionBe
  real, save :: sim_tradBe

  real, save :: sim_vacThickness
  real, save :: sim_rhoVa
  real, save :: sim_teleVa
  real, save :: sim_tionVa
  real, save :: sim_tradVa

end module Simulation_data


