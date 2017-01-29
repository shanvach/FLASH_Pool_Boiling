!!****if* source/Simulation/SimulationMain/INavierStokes/2D/LidDrivenCavity/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  Simulation_data()
!!
!! DESCRIPTION
!!
!!  Stores the local data for Simulation setup: INS-iso-turb
!!
!!***

module Simulation_data

  implicit none

#include "constants.h"

  !! *** Runtime Parameters *** !!
  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax
  logical, save :: sim_gCell
  !integer, save :: stretching_exp=1;
  integer, save :: sim_meshMe
  integer,save :: Npmax, Npmin, Ypmax
  real,save:: per=0.02; 
  integer, allocatable, save :: Nplus(:), Nminus(:)
  integer, save :: pt_Numpart
  real, save :: domainsize
end module Simulation_data
