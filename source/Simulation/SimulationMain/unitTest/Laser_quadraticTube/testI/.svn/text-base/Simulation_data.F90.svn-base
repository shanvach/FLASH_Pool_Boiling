!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testI/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data 
!!
!! DESCRIPTION
!!
!!  Stores the local data for the laser quadratic tube unit test.
!!  
!!***

Module Simulation_data

  implicit none

#include "constants.h"

  character (len = MAX_STRING_LENGTH), save :: sim_baseName

  logical, save :: sim_printBlockVariables
  logical, save :: sim_XZtypeRays

  integer, save :: sim_geometry
  integer, save :: sim_globalComm
  integer, save :: sim_globalMe
  integer, save :: sim_nFocalRays
  integer, save :: sim_refinementLevel

  real,    save :: sim_beamFrequency
  real,    save :: sim_beamWavelength
  real,    save :: sim_nc
  real,    save :: sim_nuw
  real,    save :: sim_nw
  real,    save :: sim_A
  real,    save :: sim_powerDecayFactor
  real,    save :: sim_R
  real,    save :: sim_symmetryTolerance
  real,    save :: sim_Tw
  real,    save :: sim_xc
  real,    save :: sim_xw
  real,    save :: sim_yfocal
  real,    save :: sim_Z
  real,    save :: sim_zc
  real,    save :: sim_zw

  integer, parameter :: GRID_3DCARTESIAN    =  1, &
                        GRID_2DCYLINDRICAL  =  2, &
                        GRID_2DCARTESIAN    =  3

  integer, parameter :: sim_nRaysMax            =  8
  integer, parameter :: sim_refinementLevelMax  =  6
  real,    parameter :: sim_keV2Kelvin          =  11604505.

  real,    save :: sim_rayPexit (1:sim_nRaysMax)
  real,    save :: sim_rayXexit (1:sim_nRaysMax)
  real,    save :: sim_rayZexit (1:sim_nRaysMax)

  real,    save :: sim_focalPercentErrorX   (1:sim_refinementLevelMax)
  real,    save :: sim_focalPercentErrorXZ  (1:sim_refinementLevelMax)
  real,    save :: sim_powerPercentErrorX   (1:sim_refinementLevelMax)
  real,    save :: sim_powerPercentErrorXZ  (1:sim_refinementLevelMax)

  real,    save :: sim_rayPexitPercentError (1:sim_nRaysMax)
  real,    save :: sim_rayXexitPercentError (1:sim_nRaysMax)
  real,    save :: sim_rayZexitPercentError (1:sim_nRaysMax)

end module Simulation_data
