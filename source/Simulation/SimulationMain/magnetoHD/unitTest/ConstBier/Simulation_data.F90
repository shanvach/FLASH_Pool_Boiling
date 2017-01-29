!!****if* source/Simulation/SimulationMain/magnetoHD/unitTest/ConstBier/Simulation_data
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
  real, save :: sim_ptot
  real, save :: sim_nele1
  real, save :: sim_nele2
  real, save :: sim_pele1
  real, save :: sim_pele2

  real, save :: sim_singleSpeciesA
  real, save :: sim_singleSpeciesZ

  real, save :: sim_xmin
  real, save :: sim_xmax
  real, save :: sim_ymin
  real, save :: sim_ymax

  !! *** Physical Constants *** !!
  real, save :: sim_avogadro
  real, save :: sim_boltzmann
  real, save :: sim_qele
  real, save :: sim_speedlt

end module Simulation_data


