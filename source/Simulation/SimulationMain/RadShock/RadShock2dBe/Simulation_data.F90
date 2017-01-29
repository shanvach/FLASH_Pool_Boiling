!!****if* source/Simulation/SimulationMain/RadShock/RadShock2dBe/Simulation_data
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

  ! Xenon Material Properties:
  real, save :: sim_rhoXe
  real, save :: sim_tionXe
  real, save :: sim_teleXe
  real, save :: sim_tradXe
  real, save :: sim_zbarXe
  real, save :: sim_abarXe

  ! Polyimide Material Properties:
  real, save :: sim_rhoCh
  real, save :: sim_tionCh
  real, save :: sim_teleCh
  real, save :: sim_tradCh
  real, save :: sim_zbarCh
  real, save :: sim_abarCh

  ! Vacuum Material Properties:
  real, save :: sim_rhoVa
  real, save :: sim_tionVa
  real, save :: sim_teleVa
  real, save :: sim_tradVa
  real, save :: sim_zbarVa
  real, save :: sim_abarVa

  ! Beryllium Material Properties:
  real, save :: sim_rhoBe
  real, save :: sim_tionBe
  real, save :: sim_teleBe
  real, save :: sim_tradBe
  real, save :: sim_zbarBe
  real, save :: sim_abarBe
  real, save :: sim_gradSize

  ! ************************
  ! * Geometric Properties *
  ! ************************

  ! Plastic tube, inner radius
  real, save :: sim_tubeRadius
  
  ! Plastic tube, wall thickness
  real, save :: sim_tubeThickness

  ! Beryllium slab thickness:
  real, save :: sim_slabThickness

  ! Beryllium slab radius:
  real, save :: sim_slabRadius

  ! Thickness of the Vacuum region between z = 0 and the slab:
  real, save :: sim_vacThickness

end module Simulation_data


