!!****if* source/Simulation/SimulationMain/LaserSlabLayers/Simulation_data
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

#include "constants.h"

  !! *** Runtime Parameters *** !!  
  real, save :: sim_targetRadius
  real, save :: sim_coatHeight
  real, save :: sim_targetHeight
  real, save :: sim_vacuumHeight

  real,    save :: sim_rhoTarg  
  real,    save :: sim_teleTarg 
  real,    save :: sim_tionTarg 
  real,    save :: sim_tradTarg 
  real,    save :: sim_zminTarg
  integer, save :: sim_eosTarg

  real,    save :: sim_rhoCoat
  real,    save :: sim_teleCoat
  real,    save :: sim_tionCoat
  real,    save :: sim_tradCoat	
  real,    save :: sim_zminCoat
  integer, save :: sim_eosCoat

  real,    save :: sim_rhoCham  
  real,    save :: sim_teleCham 
  real,    save :: sim_tionCham 
  real,    save :: sim_tradCham 
  integer, save :: sim_eosCham  

  real, save :: sim_smallX

end module Simulation_data


