!!****if* source/Simulation/SimulationMain/unitTest/Burn/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!    use Simulation_data
!!
!! DESCRIPTION
!!
!!  The data module for storing simulation initialization data
!!
!! ARGUMENTS
!!
!! PARAMETERS
!!
!!     smallx       smallest allowable abundance
!!     xmin         boundary of domain
!!     xmax         boundary of domain
!!     ymin         boundary of domain
!!     ymax         boundary of domain
!!     zmin         boundary of domain
!!     zmax         boundary of domain
!!     t_min        lowest temperature in domain
!!     t_max        highest temperature in domain
!!     rho_min      lowest density in domain
!!     rho_max      highest density in domain
!!     compa        name of material at left of domain
!!     compb        name of material at right of domain
!!
!!  NOTES
!!     These internal but shared variables are also defined
!!     compIndexA     species index of material at left of domain (within SPECIES_BEGIN to SPECIES_END)
!!     compIndexB     species index of material at right of domain (within SPECIES_BEGIN to SPECIES_END)
!!
!!***

Module Simulation_data
  
  implicit none

  real, save     :: sim_rhoMin, sim_rhoMax
  real, save     :: sim_tempMin,   sim_tempMax
  real, save     :: sim_imax, sim_imin, sim_jmax, sim_jmin, sim_kmax, sim_kmin
  character (len=4) :: sim_compA, sim_compB
  integer, save   :: sim_compIndexA, sim_compIndexB
  real, save     :: sim_smallx

end Module Simulation_data


