module Simulation_data

  ! Parameters of the problem.

  character(len=80) :: sim_profile1, sim_profile2
  real, save    :: sim_refiningRadius, sim_deRefiningRadius
  real, save    :: sim_refinementDensityCutoff
  integer, save :: sim_subZones, sim_lrefineMin, sim_lrefineMax, sim_ptdirn
  real, save    :: sim_smalle, sim_smallP
  real, save    :: sim_plasmaBeta
  real, save    :: sim_smallT, sim_pi
  real, save    :: sim_d
  real, save    :: sim_xMax, sim_xMin, sim_yMax, sim_yMin
  real, save    :: sim_zMax, sim_zMin, sim_smallX, sim_smlRho
  logical, save :: sim_testSingleGalaxy, sim_testAtmosphere
  logical, save :: sim_restart, sim_forceHydroLimit
  logical, save :: sim_pressureNormalize
  logical, save :: sim_useMeKaLCooling
  real, save    :: sim_Bmag, sim_lMax, sim_lMin
  
  ! Other variables
  
  integer, save :: lbox, niq

  real, save    :: nsubinv, nsubvolinv
  real, save    :: sim_xCtr, sim_yCtr, sim_zCtr
  real, save    :: sim_Newton, sim_vInit, sim_zSol, sim_Bavg
  real, save    :: sim_vrCtr, sim_xExtent, sim_yExtent, sim_zExtent
  real, save    :: sim_arCtr, sim_rCtr
  real, save    :: sim_oarCtr, sim_densAmbient, sim_presAmbient
  integer, save :: sim_nBzones
  real, save :: sim_Bxmin, sim_Bymin, sim_Bzmin, sim_Bdx, sim_Bdy, sim_Bdz
  real, save :: sim_BLx, sim_BLy, sim_BLz
  real, allocatable, dimension(:), save :: sim_Bxcoord, sim_Bycoord, &
       sim_Bzcoord
  real, save :: N_a, k_B, m_e, mueinv, sim_velDx, sim_velDy, sim_velDz

  logical, save :: sim_killdivb, sim_cleaningDivB

  real, save, dimension(:,:,:), allocatable :: sim_Bx
  real, save, dimension(:,:,:), allocatable :: sim_By
  real, save, dimension(:,:,:), allocatable :: sim_Bz

  integer, save :: sim_meshMe, numPoints1, numPoints2

  real, allocatable, save, dimension(:) :: r1, dens1, pres1, gpot1, grav1
  real, allocatable, save, dimension(:) :: r2, dens2, pres2, gpot2, grav2

end module Simulation_data
