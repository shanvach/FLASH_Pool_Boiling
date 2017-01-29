!!****if* source/Simulation/SimulationMain/magnetoHD/Riemann2D/Simulation_data
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
!!  Stores the local data for Simulation setup
!!
!!***

module Simulation_data

  implicit none

#include "constants.h"

  !! *** Runtime Parameters *** !!    
  real, save   :: sim_gamma
  real, save   :: sim_dens1,sim_pres1,sim_velx1,sim_vely1,sim_velz1,sim_magx1,sim_magy1,sim_magz1
  real, save   :: sim_dens2,sim_pres2,sim_velx2,sim_vely2,sim_velz2,sim_magx2,sim_magy2,sim_magz2
  real, save   :: sim_dens3,sim_pres3,sim_velx3,sim_vely3,sim_velz3,sim_magx3,sim_magy3,sim_magz3
  real, save   :: sim_dens4,sim_pres4,sim_velx4,sim_vely4,sim_velz4,sim_magx4,sim_magy4,sim_magz4

  !! Simulation variables
  real, save    :: sim_smallx, sim_smallP
  logical, save :: sim_gCell, sim_killdivb

  integer, save :: sim_meshMe
end module Simulation_data
