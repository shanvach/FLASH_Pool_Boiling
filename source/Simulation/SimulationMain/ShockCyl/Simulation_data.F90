!!****if* source/Simulation/SimulationMain/ShockCyl/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!  module Simulation_data
!!
!! DESCRIPTION
!!  Store the simulation data for Shock Cylinder setups
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
  include 'sim_rzDatafileSizes.fh'


  integer, save ::                          sim_meshMe
  
  integer, save ::                          numRefineVars
  integer, dimension(MAXREFVARS) , save ::  refine_var
  real, dimension(MAXREFVARS), save ::      dref
  
  real, save    :: mach, ref_rect_x, ref_rect_y, xctr, yctr
  real, save    :: xmin, xmax, ymin, ymax, zmin, zmax


  real, save    :: rho_amb, p_amb, vx_amb
  real, save    :: mw_air, mw_sf6, gamma_air, gamma_sf6
  real, save    :: rho_ps, p_ps, vx_ps
  real, save    :: sim_xShock
  real, save    :: zctr
  real, save    :: sim_smallSF6  ! minimum value for SF6 concentration
  logical, save :: sim_useRawData, use_rz_sim_data, sim_useRadialFit
  real, save    :: maxconc,  c_initial
  real, save    :: ximgmax, yimgmax
  real,       dimension(608,468) :: vconc 
  real, save, dimension(608)     :: sizex 
  real, save, dimension(468)     :: sizey
  real, save    :: sim_radialFitRadius,  sim_rawPixelSize, vz_fact
  logical, save :: rz_3d_use_sym = .false.
  real, save    :: rz_pert_zlen, rz_pert_amp, pi
  real, save    :: rz_rmax, rz_zmax, rz_zplane
  real, save    :: dr_rz, dz_rz, dri, dzi
  integer, save :: rz_subintNX, rz_subintNY, rz_subintNZ

!! EOS stuff
  real, save, dimension(EOS_NUM) :: eos_arr
!  Was set up as this but we write to MASK below the lower bound
!!  logical,save , dimension(EOS_VARS+1:EOS_NUM) :: mask
  integer, save :: mode
  integer, save :: nr_c, nz_c, nr_e, nz_e


!! For lattice initialization, only with particles
#ifdef FLASH_PARTICLES
  real, save    :: rz_xpmin, rz_xpmax, rz_ypmin, rz_ypmax, rz_zpmin, rz_zpmax
  integer, save :: rz_nxp, rz_nyp, rz_nzp
#endif

  real, save    :: r_c(nr_c_max), r_e(0:nr_c_max)
  real, save    :: z_c(nz_c_max), z_e(0:nz_c_max)
  character(len=MAX_STRING_LENGTH), save  :: sf6_file, &
       press_file, rvel_file, zvel_file

  logical, save :: gcell=.true.

  real, save    :: x_sf6_cc(  nr_c_max,   nz_c_max),  &
                   x_sf6_r1(  nr_c_max,   nz_c_max),  &
                   x_sf6_r2(  nr_c_max,   nz_c_max),  &
                   x_sf6_z1(  nr_c_max,   nz_c_max),  &
                   x_sf6_z2(  nr_c_max,   nz_c_max)

  real, save    ::  rvel_ec(0:nr_c_max,   nz_c_max),  &
                    rvel_r1(0:nr_c_max,   nz_c_max),  &
                    rvel_r2(0:nr_c_max,   nz_c_max),  &
                    rvel_z1(0:nr_c_max,   nz_c_max),  &
                    rvel_z2(0:nr_c_max,   nz_c_max)

  real, save    ::  zvel_ce(  nr_c_max, 0:nz_c_max),  &
                    zvel_r1(  nr_c_max, 0:nz_c_max),  &
                    zvel_r2(  nr_c_max, 0:nz_c_max),  &
                    zvel_z1(  nr_c_max, 0:nz_c_max),  &
                    zvel_z2(  nr_c_max, 0:nz_c_max)


  integer, save :: species_sf6
  integer, save :: species_air

  integer, save :: ir, iz, ir_c, iz_c, ir_e, iz_e

  real, save :: vmax, sim_rawMinX, sim_rawMinY
  real, save :: vz_sf6, xn_sf6
  integer, save    :: sim_rawNumPixelsX, sim_rawNumPixelsY

! defined only for paramesh 
#ifdef FLASH_GRID_PARAMESH
  integer, save :: nrefs
#endif

end module Simulation_data
