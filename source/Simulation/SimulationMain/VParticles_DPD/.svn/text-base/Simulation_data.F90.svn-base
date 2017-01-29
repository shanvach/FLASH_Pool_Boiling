!!****if* source/Simulation/SimulationMain/VParticles_DPD/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!  module Simulation_data()
!!
!! DESCRIPTION
!!
!!  Store the simulation data for unitTesting of Particles
!!   
!! ARGUMENTS
!!
!! PARAMETERS   
!!      Described in the Config file
!!
!! NOTES
!!
!!  No arguments.  All data passed by "use Simulation_data"
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"
#include "Flash.h"
  
  real, save    :: sim_rho_amb, sim_p_amb, sim_vx_amb, sim_vx_multiplier
  real, save    :: sim_seed, sim_vx_pert, sim_vy_pert, sim_vz_pert
  real, dimension(MDIM), save :: sim_initPos, sim_deltaMove
  integer, save :: sim_meshMe
  real, save    :: ks,kd,rc,lo,lambda,domainsize(MDIM)
  real, save    :: sig,sig_sq,kt
  integer,save  ::  pt_NumBTypes
  integer,save  ::  pt_NumPart, pt_NumBodies,pt_BodyTypes
  logical,save  ::  Firstcall=.true.
  integer,save,allocatable :: beadsPBody(:)
  Type Conn
   integer :: numLinks   
   integer, allocatable, dimension(:,:) ::  Links
  end Type Conn

  Type(Conn), allocatable, dimension(:) :: Connect 
  integer,parameter :: MAXBTYPES=10;
  real,save :: Aij(MAXBTYPES,MAXBTYPES)
  real,save :: sqrt_dt
  integer,parameter:: MAXCONN=10
  real,save :: eps ,Aslab,local_Aslab,local_surftension
end module Simulation_data
