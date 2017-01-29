!!****if* source/Simulation/SimulationMain/PoisTest_particleBasedRefine/Particles_initPositions
!!
!! NAME
!!    Particles_initPositions
!!
!! SYNOPSIS
!!
!!    Particles_initPositions( logical, INTENT(inout) :: success)
!!
!! DESCRIPTION
!!    Initialize particle locations.  This version is designed for the 
!!    particle based poistest problem.  The particle locations are
!!    initialized according to the density given by a 
!!    Gaussian distribution function.
!!
!! ARGUMENTS
!!
!!    success : boolean indicating whether particles positions were 
!!              successfully initialized. This is not really relevant
!!              for this version of the routine
!!
!!
!!***

subroutine Particles_initPositions (success,updateRefine)

  use Particles_data, ONLY:  pt_numLocal, particles, pt_maxPerProc, &
       pt_numAtOnce,pt_posInitialized

  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkBoundBox, &
       Grid_getBlkPtr, Grid_releaseBlkPtr
  use pt_initFromFileInterface, ONLY : pt_initNumToGet, pt_initNextNParticles, pt_initIfInBlock

  use Simulation_interface, ONLY : Simulation_initBlock


  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  logical, intent(INOUT) :: success
  logical, INTENT(OUT) :: updateRefine

  integer       :: i, j, k, n, m, ierr
  integer       :: pb,pe
  integer       :: blockID
  logical       :: notDone,localSuccess
  integer       :: blkCount, numToGet,numReturned
  integer,dimension(MAXBLOCKS) :: blkList
  real, dimension(LOW:HIGH,MDIM) :: boundBox
  integer, save :: numGot=0
  real,dimension(:,:,:,:),pointer :: solnData
  !----------------------------------------------------------------------


  interface
     subroutine InitializeParticles()
     end subroutine InitializeParticles
  end interface


  print *, "In Particles_initPositions, pt_posInitialized=", pt_posInitialized 

  if(pt_posInitialized) then
     success = .true.
     return
  end if
  pt_numLocal=0

  !We must determine the initial number of top level blocks.
  call Grid_getListOfBlocks(LEAF,blkList,blkCount)


  !Initialize the particles based on the density distribution in DENS_VAR.
  call InitializeParticles()


  !Now that we have initialized the particles, we are no longer interested in 
  !DENS_VAR.  Therefore, store a constant small density in the DENS_VAR grid 
  !element, so that the gravitational potential is not affected.
  do i = 1, blkCount
     call Grid_getBlkPtr(blkList(i), solnData)
     solnData(DENS_VAR,:,:,:) = 0.000001
     call Grid_releaseBlkPtr(blkList(i), solnData)
  end do


  success = .true.
  pt_posInitialized = success
  return

end subroutine Particles_initPositions
