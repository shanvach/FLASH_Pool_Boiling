!!****if* source/Simulation/SimulationMain/VParticles_DPD/Particles_initPositions
!!
!! NAME
!!    Particles_initPositions
!!
!! SYNOPSIS
!!
!!    Particles_initPositions( logical, INTENT(inout) :: partPosInitialized,
!!                             logical, INTENT(out) :: updateRefine)
!!
!!
!! DESCRIPTION
!!
!!    Initialize particle locations. This routine calls pt_initPositions
!!    which a Particles unit's local routine to initialize the positions
!!    on leaf blocks. The routine also creates tags for all the particles
!!    This routine will initialize based on Lattice or with Density 
!!    distribution  depending upon which of the two is selected. 
!!
!! ARGUMENTS
!!
!!  partPosInitialized : boolean indicating whether particles positions were 
!!            successfully initialized. This is not really relevant
!!            for this version of the routine
!!
!! updateRefine : is true if the routine wished to retain the already
!!                initialized particles instead of reinitializing them
!!                as the grid refine.
!!
!!***

!#define DEBUG_PARTICLES

subroutine Particles_initPositions (partPosInitialized,updateRefine)


  use Grid_interface, ONLY : Grid_getListOfBlocks,Grid_getBlkBoundBox
  use Driver_interface, ONLY : Driver_abortFlash
  use pt_interface, ONLY : pt_initPositions,pt_createTag
  use Particles_data, ONLY : pt_posInitialized,pt_numLocal,useParticles,&
       pt_typeInfo, particles, pt_meshNumProcs, pt_posInitialized, pt_meshMe
  use Simulation_data,ONLY : sim_initPos
  use pt_interface, ONLY :  pt_initLocal
 
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Particles.h"

  logical, INTENT(INOUT) :: partPosInitialized
  logical, INTENT(OUT) :: updateRefine

  integer       :: i, j, k, b
  integer       :: p
  integer       :: numLocalThisType, numNewLocalThisType, numLocalPreviousTypes
  integer       :: numPreviousLocal
  integer       :: blockID
  logical       :: IsInBlock
  real          :: xpos, ypos, zpos, bxl, byl, bzl, bxu, byu, bzu
  real          :: xvel, yvel, zvel

! NOTE dxParticle is particle spacing, not grid spacing
  real, dimension(MDIM) :: dxParticle = 0.0
  real, dimension(2,MDIM):: boundBox
  integer :: blkCount
  integer,dimension(MAXBLOCKS) :: blkList
!----------------------------------------------------------------------

  if(.not.useParticles) then
     partPosInitialized = .true.
  end if
  if(partPosInitialized) return

  !CD: We need to move the call to pt_initLocal to this level. 
  !Otherwise, if it is in pt_initPositions, we get a deadlock
  !when the number of blocks are not the same on all processors.
  call pt_initLocal()

  updateRefine = .false.
  pt_numLocal = 0

  !CD: The code in pt_initPositions should be moved into this subroutine.
  call pt_initPositions(-1,partPosInitialized)
  if (.not.partPosInitialized) then
     call Driver_abortFlash("pt_initPositions did not initialize particles")
  end if
  pt_posInitialized = partPosInitialized

!!$   write(*,*) 'pt_numlocal',pt_numlocal
!!$  do i= 1,pt_numlocal
!!$     write(*,*) i,particles(BLK_PART_PROP,i)
!!$  end do
  call pt_createTag()

end subroutine Particles_initPositions
