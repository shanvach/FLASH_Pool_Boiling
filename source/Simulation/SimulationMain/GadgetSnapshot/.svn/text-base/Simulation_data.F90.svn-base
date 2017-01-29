module Simulation_data

  ! Runtime parameters
 
  character(len=80), save :: sim_path, sim_basename
  integer, save :: sim_numParticles
  integer, save :: sim_numFiles, sim_snapshotNumber
  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax
  real, save    :: sim_zMin, sim_zMax

  ! Other variables

  logical, save :: sim_useCosmo = .false.
  real, save    :: sim_unitLength, sim_unitTime, sim_unitVelocity
  real, save    :: sim_unitMass, sim_unitEnergy, sim_unitPressure
  real, save    :: sim_hubble, sim_omegaMatter, sim_omegaLambda
  real, save    :: sim_boxSize, sim_time, sim_redshift, sim_useCosmo
  real, save, dimension(:), allocatable :: posx, posy, posz
  real, save, dimension(:), allocatable :: velx, vely, velz
  real, save, dimension(:), allocatable :: mass, type, tags

  integer, save :: sim_meshMe
end module Simulation_data
