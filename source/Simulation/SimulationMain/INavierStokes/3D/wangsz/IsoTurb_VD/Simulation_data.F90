!!****if* source/Simulation/SimulationMain/INavierStokes/3D/kientzler_bubbleFS_caseA/Simulation_data
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
  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
!  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax
  logical, save :: sim_gCell

end module Simulation_data
