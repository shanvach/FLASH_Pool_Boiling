!!****if* source/Simulation/SimulationMain/AdvectMassScalars/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!  use Simulation_data
!!
!! DESCRIPTION
!!  Store the simulation data for advection of mass scalars setups
!!   
!! ARGUMENTS
!!  None.  All data passed by "use Simulation_data"
!!
!! PARAMETERS   
!!      Described in the Config file
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"
#include "Eos.h"
#include "Flash.h"

  real, save :: sim_xcos, sim_ycos, sim_zcos

  real, save :: sim_smallp, sim_smallx

  real, save :: sim_gamma

  real, save :: sim_rhoin, sim_rhoout, sim_msin, sim_msout, sim_pressure, sim_velocity, &
       sim_width, sim_phase, sim_xangle, sim_yangle, sim_zangle, sim_posn

  integer, save :: sim_pulse_fctn, sim_pulse_fctn_ms1, sim_pulse_fctn_ms2, sim_pulse_fctn_ms3, &
                   sim_pulse_fctn_ms4, sim_pulse_fctn_ms5

  integer, save :: sim_ims1, sim_ims2, sim_ims3, sim_ims4, sim_ims5
  logical, save :: sim_planar
  real, save    :: sim_xmin, sim_xmax, sim_ymin, sim_ymax, sim_zmin, sim_zmax, &
                   sim_pi, sim_twoPi, sim_multid, sim_xcent, sim_ycent, sim_zcent

  integer, save :: sim_meshMe
end module Simulation_data
