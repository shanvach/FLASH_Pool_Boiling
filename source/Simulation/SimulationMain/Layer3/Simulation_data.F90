!!****if* source/Simulation/SimulationMain/Layer3/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!  use Simulation_data
!!
!! DESCRIPTION
!!  Store the simulation data for Three layer setup
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
      
  integer, parameter ::  N1D = 8192

  real, save, dimension(N1D,NUNK_VARS) :: sim_1dModel
  real, save, dimenstion(N1D) :: sim_xzn1d

  real, save :: sim_pi

  real, save :: sim_xmin, sim_xmax, sim_ymin, sim_ymax
  real, save :: sim_small
  real, save :: sim_gamcu, sim_gamcf, sim_gamch

      integer, PARAMETER :: & 
     &     sim_dens = 1,                      & 
     &     sim_velx = 2,                                       & 
     &     sim_vely = 3,                                       & 
     &     sim_velz = 4,                                       & 
     &     sim_pres = 5,                                       & 
     &     sim_ener = 6,                                       & 
     &     sim_temp = 7,                                       & 
     &     sim_gamc = 8,                                       & 
     &     sim_game = 9,                                       & 
     &     sim_enuc = 10,                                      & 
     &     sim_gpot = 11, & 
     &     sim_speciesBegin = 12                

      
  integer, save :: sim_meshMe
end module Simulation_data


