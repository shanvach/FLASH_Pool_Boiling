!!****if* source/Simulation/SimulationMain/Nuc2Grid/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!!
!! DESCRIPTION
!!
!!  Initializes all the parameters needed for the Huang & Greengard (poistest)
!!  problem
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!  sim_smlRho : 
!!  sim_gam : 
!!***

subroutine Simulation_init()
  
  use Simulation_data, ONLY : sim_ptMass, sim_densityThreshold,sim_nucFileNames,&
       sim_smlRho,sim_doWAvg,sim_doGP,sim_doWeight,sim_doConvolve,sim_doInterpExtrap,&
       sim_doLowerBounds,sim_doEos, &
       sim_smallE, sim_smallT, sim_meshMe, sim_meshNumProcs, &
       sim_convoSmearWidI, sim_convoSmearWidJ, sim_convoSmearWidK, &
       sim_convoSmearShapeI, sim_convoSmearShapeJ, sim_convoSmearShapeK, &
       sim_ptInNdim, sim_ptInGeometryStr, sim_geometryStr, sim_ptInGeometry, sim_geometry, &
       sim_MaxParticleFiles, sim_ptNumPartFiles, &
       sim_doFixupAbundances, sim_abundanceFixupMaxDens, &
       sim_restart, sim_unkCellWeight, sim_useTrajValues
  use Particles_data, ONLY: pt_numParticlesWanted

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, RuntimeParameters_mapStrToInt
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype, Driver_getNumProcs
  use Logfile_interface, ONLY : Logfile_stamp
  implicit none
#include "Flash.h"
#include "constants.h"


  integer :: i
  character(len=MAX_STRING_LENGTH) :: rpName !DEV: Do we need SAVE here? - assume not. KW

  call Logfile_stamp("initializing for NucOToRT data converter",  &
       "[Simulation_init]")
  call RuntimeParameters_get("restart",sim_restart)
  call RuntimeParameters_get("sim_ptMass",sim_ptMass)
  call RuntimeParameters_get("sim_densityThreshold",sim_densityThreshold)
!  sim_densityThreshold = 0.0
  call RuntimeParameters_get("sim_smlRho",sim_smlRho)
  
  call RuntimeParameters_get("sim_ptNumPartFiles",sim_ptNumPartFiles)
  if (sim_ptNumPartFiles > sim_MaxParticleFiles) then
     call Driver_abortFlash("Too many particle files, increase sim_MaxParticleFiles!")
  end if
  do i = 1,sim_ptNumPartFiles
     call concatStringWithInt("sim_nucFileName_",i,rpName)
     call RuntimeParameters_get(rpName,sim_nucFileNames(i))
  end do
  call RuntimeParameters_get("pt_numParticlesWanted", pt_numParticlesWanted)

  call RuntimeParameters_get("doWeight",sim_doWeight)
  call RuntimeParameters_get("doWAvg",sim_doWAvg)
  call RuntimeParameters_get("doGP",sim_doGP)
  call RuntimeParameters_get("doConvolve",sim_doConvolve)
  call RuntimeParameters_get("doInterpExtrap",sim_doInterpExtrap)
  call RuntimeParameters_get("doLowerBounds",sim_doLowerBounds)
  call RuntimeParameters_get("smalle",sim_smallE)
  call RuntimeParameters_get("smallt",sim_smallT)
  call RuntimeParameters_get("doEos",sim_doEos)
  call RuntimeParameters_get("doFixupAbundances",sim_doFixupAbundances)
  call RuntimeParameters_get("sim_abundanceFixupMaxDens",sim_abundanceFixupMaxDens)

  call RuntimeParameters_get("convoSmearWidI",sim_convoSmearWidI)
  call RuntimeParameters_get("convoSmearWidJ",sim_convoSmearWidJ)
  call RuntimeParameters_get("convoSmearWidK",sim_convoSmearWidK)
  call RuntimeParameters_get("convoSmearShapeI",sim_convoSmearShapeI)
  call RuntimeParameters_get("convoSmearShapeJ",sim_convoSmearShapeJ)
  call RuntimeParameters_get("convoSmearShapeK",sim_convoSmearShapeK)

  call RuntimeParameters_get("particlesInputNdim",sim_ptInNdim)
  call RuntimeParameters_get("particlesInputGeometry",sim_ptInGeometryStr)
  call RuntimeParameters_mapStrToInt(sim_ptInGeometryStr, sim_ptInGeometry)
  call RuntimeParameters_get("geometry",sim_geometryStr) !geometry for Grid
  call RuntimeParameters_mapStrToInt(sim_geometryStr, sim_geometry)

  call Driver_getMype(MESH_COMM, sim_meshMe)
  call Driver_getNumProcs(MESH_COMM, sim_meshNumProcs)

  if (sim_restart) then
     ! dens,eint,gamc*,game*,temp,velx,vely,velz; *=irrelevant
     sim_useTrajValues(PROP_VARS_BEGIN:PROP_VARS_END) = .FALSE.
     sim_unkCellWeight(PROP_VARS_BEGIN:PROP_VARS_END) = 1.0
     sim_unkCellWeight(NUMP_VAR) = 0.0 !but special in Driver_evolveFlash anyway
     sim_unkCellWeight(NUP0_VAR) = 0.0 !but special in Driver_evolveFlash anyway
     sim_unkCellWeight(NUP1_VAR) = 0.0 !but special in Driver_evolveFlash anyway
     sim_unkCellWeight(PDEN_VAR) = 0.0 !but exempted in Driver_evolveFlash anyway
     sim_unkCellWeight(ENTR_VAR) = 0.0 !none there
     sim_useTrajValues(PRES_VAR) = .false. !none there
     sim_unkCellWeight(PRES_VAR) = 1.0 !but may get recomputed

     sim_useTrajValues(SPECIES_BEGIN:SPECIES_END) = .true.
     sim_unkCellWeight(SPECIES_BEGIN:SPECIES_END) = 0.0
!!$     sim_unkCellWeight(C_SPEC) = 1.0e-21
!!$     sim_unkCellWeight(O_SPEC) = 1.0e-21

     sim_useTrajValues(MASS_SCALARS_BEGIN:MASS_SCALARS_END) = .true.

     sim_useTrajValues(RPV1_MSCALAR) = .false.
     sim_useTrajValues(RPV2_MSCALAR) = .false.
     sim_useTrajValues(RPV3_MSCALAR) = .false.

     sim_unkCellWeight(MASS_SCALARS_BEGIN:MASS_SCALARS_END) = 1.0
     sim_useTrajValues(SUMX_MSCALAR) = .false. !none there
     sim_unkCellWeight(SUMY_MSCALAR) = 0.0
     sim_unkCellWeight(YE_MSCALAR) = 0.0
     sim_unkCellWeight(NI56_MSCALAR) = 0.0

  else
     sim_useTrajValues(1:NUNK_VARS) = .true.
     sim_unkCellWeight(1:NUNK_VARS) = 0.0
  end if
  where(sim_unkCellWeight .NE. 0.0 .AND. .NOT. sim_useTrajValues)
     sim_unkCellWeight = 1.0
  end where

!!  sim_useAnyUnkValues = (ANY(sim_unkCellWeight(:).NE.0.0))

  call sim_initOutputGrid()

end subroutine Simulation_init
