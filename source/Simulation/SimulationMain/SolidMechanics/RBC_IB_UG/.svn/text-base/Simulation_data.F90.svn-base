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

  integer, save :: sim_meshMe
  real, save :: rc
  integer,save :: Npmax, Npmin, Ypmax
  real,save:: per=0.02;
  integer, allocatable, save :: Nplus(:), Nminus(:)

end module Simulation_data
