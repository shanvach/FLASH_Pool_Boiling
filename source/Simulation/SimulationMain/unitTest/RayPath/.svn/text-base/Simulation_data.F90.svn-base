!!****if* source/Simulation/SimulationMain/unitTest/RayPath/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data 
!!
!!  DESCRIPTION
!!
!!  Stores the local data for the unit test of path of a ray in presence 
!!  of refraction 
!!  
!! PARAMETERS
!!
!!  sim_ilBnd         The lower IAXIS bound of the refracting slab
!!  sim_iubnd         The upper IAXIS bound 
!!  sim_jlBnd         The lower JAXIS bound of the refracting slab
!!  sim_jubnd         The upper JAXIS bound 
!!  sim_klBnd         The lower KAXIS bound of the refracting slab
!!  sim_kubnd         The upper KAXIS bound 
!!  sim_refract       Refractive index of the slab
!!  sim_numRay        The number of rays 
!!  sim_refractType   the type of refraction (constant, or with constant
!!                    gradient)
!!  sim_rayIncidence  array containing the initial entry point and incidence
!!                    angles of each ray.
!!  sim_fileRay       filename of the file containing the incidence info
!!
!!***

module Simulation_data

  implicit none

  !! *** Runtime Parameters *** !!

  real, save, dimension(LOW:HIGH,MDIM) :: sim_slabBndBox,sim_domainBndBox
  real, save :: sim_ilBnd, sim_iuBnd, sim_jlBnd, sim_juBnd
  real, save    :: sim_klBnd, sim_kuBnd
  real, save    :: sim_refract
  integer, save :: sim_numRay
  real, save, allocatable(:,:) :: sim_rayIncidence, sim_rayFinal
  charater(len=40),save :: sim_fileRay, sim_refractType

integer, save :: sim_meshMe
end module Simulation_data
