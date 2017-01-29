!!****if* source/Simulation/SimulationMain/RadShock/RadShock2dFull/Simulation_data
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
  real, save :: sim_smallX

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

  ! Gold Material Properties:
  real, save :: sim_rhoAu
  real, save :: sim_tionAu
  real, save :: sim_teleAu
  real, save :: sim_tradAu
  real, save :: sim_zbarAu
  real, save :: sim_abarAu

  ! Acrylic Material Properties:
  real, save :: sim_rhoAc
  real, save :: sim_tionAc
  real, save :: sim_teleAc
  real, save :: sim_tradAc
  real, save :: sim_zbarAc
  real, save :: sim_abarAc

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

  ! Gold washer thickness
  real, save :: sim_goldThickness

  ! Gold washer outer radius
  real, save :: sim_goldRadius

  ! Acrylic collar thickness
  real, save :: sim_acrylicThickness

  ! Acrylic outer radius
  real, save :: sim_acrylicRadius

  ! Length of the "window" in the tube
  real, save :: sim_windowThickness

  ! ***********************************
  ! * Localized Refinement Properties *
  ! ***********************************

  ! Base lrefine_max:
  integer, save :: sim_lrefmaxBase

  ! Parameters controlling lrefmax in Beryllium region:
  integer, save :: sim_lrefmaxBe
  real   , save :: sim_belrXMin
  real   , save :: sim_belrXMax
  real   , save :: sim_belrYMin
  real   , save :: sim_belrYMax
  
  ! Parmaeters controlling lrefmax in Polyimide region:
  integer, save :: sim_lrefmaxPoly
  real   , save :: sim_polylrXMin
  real   , save :: sim_polylrXMax
  real   , save :: sim_polylrYMin
  real   , save :: sim_polylrYMax

end module Simulation_data


