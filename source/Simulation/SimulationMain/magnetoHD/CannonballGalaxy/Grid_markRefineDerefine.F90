!!****if* source/Simulation/SimulationMain/magnetoHD/CannonballGalaxy/Grid_markRefineDerefine
!!
!! NAME
!!  Grid_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  Grid_markRefineDerefine()
!!  
!! DESCRIPTION 
!!  Mark blocks for refinement or derefinement
!!  This routine is used with AMR only where individual 
!!  blocks are marked for refinement or derefinement based upon
!!  some refinement criterion. The Uniform Grid does not need
!!  this routine, and uses the stub.
!!
!! ARGUMENTS
!! 
!! NOTES
!!
!! Every unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. For Grid unit these variables begin with "gr_"
!! like, gr_meshMe or gr_eosMode, and are stored in fortran
!! module Grid_data (in file Grid_data.F90). The other variables
!! are local to the specific routines and do not have the prefix "gr_"
!!
!!
!!***

subroutine Grid_markRefineDerefine()
  
  use Simulation_data, ONLY : sim_refiningRadius, sim_deRefiningRadius, &
       sim_refinementDensityCutoff, sim_xCtr, sim_yCtr, sim_zCtr

  use Grid_data, ONLY : gr_refine_cutoff, gr_derefine_cutoff,&
                        gr_refine_filter,&
                        gr_numRefineVars,gr_refine_var,&
                        gr_refineOnParticleCount, &
                        gr_minParticlesPerBlk, gr_maxParticlesPerBlk
  use tree, ONLY : newchild, refine, derefine, stay, lrefine_min, &
       lrefine_max, lrefine
  use Grid_interface, ONLY : Grid_fillGuardCells, &
    Grid_getListOfBlocks, Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_markBlkRefine, Grid_getBlkIndexLimits, Grid_markRefineSpecialized, &
    Grid_markBlkDerefine, Grid_getCellCoords, Grid_mapParticlesToMesh, &
    Grid_getBlkPhysicalSize, Grid_getBlkCenterCoords

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  real :: ref_cut,deref_cut,ref_filter
  integer       :: l,i,iref, ierr
  logical :: doEos=.true.
  integer,parameter :: maskSize = NUNK_VARS+NDIM*NFACE_VARS
  integer :: lb, blockCount, blockID
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MAXBLOCKS) :: blockList
  logical,dimension(maskSize) :: gcMask

  real, dimension(:,:,:,:), pointer :: solnData

  real, dimension(MDIM) :: bsize

  real :: radius
  real :: localMinEnergy(2), globalMinEnergy(2), minPosn(3), specs(4), minEnergy
  real, dimension(:), allocatable :: xCoord, yCoord, zCoord 
  integer :: minPE, minIndex(3), minBlock, size(3)

  real, dimension(MDIM) :: bcoords

  integer :: oneBlkCount
  logical :: deref

  gcMask=.false.
  do i = 1,gr_numRefineVars
     iref = gr_refine_var(i)
     if (iref > 0) gcMask(iref) = .TRUE.
  end do

  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,doEos=.true.,&
       maskSize=maskSize, mask=gcMask, makeMaskConsistent=.true.)

  newchild(:) = .FALSE.
  refine(:)   = .FALSE.
  derefine(:) = .FALSE.
  stay(:)     = .FALSE.

  do l = 1,gr_numRefineVars
     iref = gr_refine_var(l)
     ref_cut = gr_refine_cutoff(l)
     deref_cut = gr_derefine_cutoff(l)
     ref_filter = gr_refine_filter(l)
     call gr_markRefineDerefine(iref,ref_cut,deref_cut,ref_filter)
  end do

#ifdef FLASH_GRID_PARAMESH2
  ! Make sure lrefine_min and lrefine_max are obeyed - KW
  if (gr_numRefineVars .LE. 0) then
     call gr_markRefineDerefine(-1, 0.0, 0.0, 0.0)
  end if
#endif

  call Grid_getListOfBlocks(LEAF, blockList, blockCount)
  
  ! Turn off second-derivative refinement in low gas-density regions

  do lb = 1, blockCount

     blockID = blockList(lb)

     call Grid_getBlkPtr(blockID, solnData)

     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)

     if (maxval(solnData(DENS_VAR, &
          blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
          blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
          blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))) < &
          sim_refinementDensityCutoff) then
        
        refine(blockID) = .false.

     endif

     call Grid_releaseBlkPtr(blockID, solnData)

  enddo

  ! Don't refine outside of a certain radius 

  do lb = 1, blockCount

     blockID = blockList(lb)

     call Grid_getBlkCenterCoords(blockID, bcoords)

     radius = sqrt((bcoords(1)-sim_xCtr)**2 + &
          (bcoords(2)-sim_yCtr)**2 + &
          (bcoords(3)-sim_zCtr)**2)

     if (radius > sim_deRefiningRadius) refine(blockID) = .false.
     
  enddo

  ! Refine a spherical volume around the center of the big cluster
  ! all the way to lrefine_max

  specs(1) = sim_xCtr
  specs(2) = sim_yCtr
  specs(3) = sim_zCtr
  specs(4) = sim_refiningRadius

  call Grid_markRefineSpecialized(INRADIUS, 4, specs, lrefine_max)

  return
end subroutine Grid_markRefineDerefine

