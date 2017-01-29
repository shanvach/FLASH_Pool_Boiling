!!****if* source/Simulation/SimulationMain/SamraiTest/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!  Store the simulation data
!!   
!!
!!***

module Simulation_data

  implicit none

  integer, save:  sim_meshMe, sim_numProcs

  !! *** Runtime Parameters *** !!
  integer, save :: sim_gridIMax,sim_gridJMax,sim_gridKMax
  integer, save :: sim_iGridSize,sim_jGridSize,sim_kGridSize

  !! *** Variables pertaining to unit test: UG *** !!
  logical, save :: sim_gCell

end module Simulation_data


