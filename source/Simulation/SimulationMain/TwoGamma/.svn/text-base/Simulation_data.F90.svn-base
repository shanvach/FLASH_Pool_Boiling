!!****if* source/Simulation/SimulationMain/TwoGamma/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!  use Simulation_data
!!
!! DESCRIPTION
!!  Store the simulation data for Two gamma setup
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

  real, save :: sim_xmin, sim_xmax, sim_ymin, sim_ymax
  real, save :: sim_small


! Parameters of the run 
  real, save :: sim_p0               ! constant pressure
  real, save :: sim_rho1             ! density of fluid 1
  real, save :: sim_rho2             ! density of fluid 2
  real, save :: sim_gam1             ! gamma of fluid 1
  real, save :: sim_gam2             ! gamma of fluid 2
  real, save :: sim_int1             ! internal energy of fluid 1
  real, save :: sim_int2             ! internal energy of fluid 2
  real, save :: sim_cvelx            ! initial velocity
  real, save :: sim_gammac1 
  real, save :: sim_gammac2          ! inital gammas of the two fluids
  real, save :: sim_xpert            ! initial perturbation
  real, save :: sim_temp1            ! temperature of fluid 1
  real, save :: sim_temp2            ! temperature of fluid 2
      
integer, save :: sim_meshMe
end module Simulation_data

