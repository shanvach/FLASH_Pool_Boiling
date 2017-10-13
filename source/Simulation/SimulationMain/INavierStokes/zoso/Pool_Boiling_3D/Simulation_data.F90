!!****if* source/Simulation/SimulationMain/INavierStokes/2D/IB_Cyl_parallelIBVP_HYPRE_VD_halfDiamBel/Simulation_data
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
  real, save    :: sim_zMin, sim_zMax
  logical, save :: sim_gCell

  integer, save :: sim_meshMe
  real, dimension(MDIM), save :: sim_initPos

  real, save :: sim_qIn1, sim_qIn2, sim_qOut1, sim_qOut2

  real, save :: sim_waveA
 
end module Simulation_data
