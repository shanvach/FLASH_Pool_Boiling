!!****if* source/Simulation/SimulationMain/Nuc2Grid/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data for the poistest problem
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!
!!   
!!
!!***

module Simulation_data
#include "Flash.h"
#include "constants.h"
  
  implicit none

  !! *** Runtime Parameters *** !!

  integer,save :: sim_meshMe, sim_meshNumProcs
  real, save :: sim_ptMass, sim_densityThreshold, sim_smlRho
  real, save :: sim_smallE, sim_smallT
  integer,parameter :: sim_MaxParticleFiles = 32
  character(len=80),save :: sim_nucFileNames(sim_MaxParticleFiles)
  integer,dimension(NPART_PROPS),save :: sim_mapToSpecies
  character(len=OUTPUT_PROP_LENGTH), save, allocatable, dimension(:) :: sim_propNames
  logical,save :: sim_doWAvg,sim_doGP,sim_doWeight
  logical, save :: sim_doConvolve,sim_doInterpExtrap,sim_doLowerBounds,sim_doEos
  logical,save :: sim_doFixupAbundances
  real,   save :: sim_abundanceFixupMaxDens
  integer,save :: sim_convoSmearShapeI,sim_convoSmearShapeJ,sim_convoSmearShapeK
  real,save :: sim_convoSmearWidI,sim_convoSmearWidJ,sim_convoSmearWidK
  integer,save :: sim_ptInNdim
  character(len=MAX_STRING_LENGTH),save :: sim_ptInGeometryStr, sim_geometryStr
  integer,save :: sim_ptInGeometry, sim_geometry
  integer,dimension(NPART_PROPS),save :: theZs, theAs, sim_pProp2A, sim_pProp2Z
  integer,save :: sim_numAbundanceProps
  character(len=24), allocatable :: theSpeciesNames(:) 

  integer,save :: sim_ptNumPartFiles
  integer :: sim_iPartFile
  real,save :: sim_unkCellWeight(NUNK_VARS)
  logical,save :: sim_restart
  logical,save :: sim_useTrajValues(NUNK_VARS)
end module Simulation_data
