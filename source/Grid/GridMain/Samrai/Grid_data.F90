!!****if* source/Grid/GridMain/Samrai/Grid_data
!!
!! NAME
!!  Grid_data
!!
!! SYNOPSIS
!!
!!  Store the grid data
!!
!!  This includes the global integer identifier for a block, the grid geometry information
!!  
!!  
!! 
!! CREATE AD:04/12/04
!!
!! 
!!   Defining data structures for storing paramesh related infomation.  
!!   including function for updating the grid information
!!
!! MODIFIED AD:05/19/04
!!   
!!***

Module Grid_data


  implicit none

#include "constants.h"
#include "Flash.h"

!! Define the block information

!  Type define

  integer, parameter :: somethingBig = 20

  real,save,dimension(3,somethingBig) :: delta
  logical, save :: fixedBlockSize = .TRUE.
  integer, save :: ilo_gc = GRID_ILO_GC
  integer, save :: ihi_gc = GRID_IHI_GC
  integer, save :: jlo_gc = GRID_JLO_GC
  integer, save :: jhi_gc = GRID_JHI_GC
  integer, save :: klo_gc = GRID_KLO_GC
  integer, save :: khi_gc = GRID_KHI_GC
  integer, save :: iguard = NGUARD
  integer, save :: jguard = NGUARD 
  integer, save :: kguard = NGUARD

  integer, save :: ilo = GRID_ILO
  integer, save :: ihi = GRID_IHI
  integer, save :: jlo = GRID_JLO
  integer, save :: jhi = GRID_JHI
  integer, save :: klo = GRID_KLO
  integer, save :: khi = GRID_KHI
  integer, parameter :: SAVED_VARS = 1

 

  !!stuff specific to samrai

  real, save :: gr_physDomainMinMax(2*MDIM)
  integer, save :: gr_indexDomainMinMax(2*MDIM)
  integer, save :: gr_iguard, gr_jguard, gr_kguard

  real, save :: gr_refineRatio
  real, save :: gr_effTolerance
  real, save :: gr_combineEff
  integer, save :: gr_indexMaxPatchSize(MDIM)
  integer, save :: gr_indexMinPatchSize(MDIM)

  integer, parameter :: gr_maxPatches = 1000000

  logical, save :: gr_justExchangedGC

  type gridBlock
     !!blockID is integer coordinates of the lower left cornor
     !! (ie the smallest point) of a block
     integer,dimension(MDIM) :: blockId
     !! atmost 2 neighbors, 2faces along
     !! each dimension, hence.
     real,dimension(3,GRID_IHI_GC) :: firstAxisCoords
     real,dimension(3,GRID_JHI_GC) :: secondAxisCoords
     real,dimension(3,GRID_KHI_GC) :: thirdAxisCoords

     integer,dimension(2,MDIM*2)  :: neighbors
     real,dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC,MDIM) :: facearea
     real,dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: cellvolume
     real,dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC,&
          NSAVED_GRID_VARS) :: saved
  end type gridBlock
  
  type(gridBlock),save,dimension(MAXBLOCKS),target :: oneBlock
  real,save, dimension(2,NYB,NZB,MAXBLOCKS) :: xarea,xdtdx,xgrav,xngrav,xfict
  real,save, dimension(NXB,2,NZB,MAXBLOCKS) :: yarea,ydtdy,ygrav,yngrav,yfict
  real,save, dimension(NXB,NYB,2,MAXBLOCKS) :: zarea,zdtdz,zgrav,zngrav,zfict
  
  integer ,save :: gr_nrefs
  logical ,save :: gr_convertToConsvdForMeshCalls,meshModified
  logical ,save :: unbiased_geometry
  logical, save :: monotone

  integer, save :: numRefineVars
  integer,allocatable,dimension(:) ,save :: refineVars
  integer ,save :: lrefine_min, lrefine_max
  real ,save :: xmin,xmax,ymin,ymax,zmin,zmax
  integer ,save :: xl_bcType
  integer ,save :: xr_bcType
  integer ,save :: yl_bcType
  integer ,save :: yr_bcType
  integer ,save :: zl_bcType
  integer ,save :: zr_bcType
  integer ,save :: geometry
  character(len=MAX_STRING_LENGTH) :: str_geometry
  integer ,save :: dominant_boundary
  integer ,save :: nBlockX, nblockY, nblockZ
  integer ,save, dimension(MAXREFVARS) :: refine_var    

  real,dimension(MAXREFVARS), save::refine_cutoff,derefine_cutoff,refine_filter
  logical, save :: msgbuffer 
  real, save :: smalle,smallrho
  integer ,dimension(MAXBLOCKS), save :: gr_blkList
end Module Grid_data





