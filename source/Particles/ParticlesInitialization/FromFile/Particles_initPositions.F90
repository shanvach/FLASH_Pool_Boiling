!!****if* source/Particles/ParticlesInitialization/FromFile/Particles_initPositions
!!
!! NAME
!!    Particles_initPositions
!!
!! SYNOPSIS
!!
!!    Particles_initPositions( logical, INTENT(out) :: partPosInitialized,
!!                             logical, INTENT(out) :: updateRefine)
!!
!! DESCRIPTION
!!    Initialize particle locations.  This version reads them in from
!!    a file. There is some error control built in here in case the
!!    particles are too closely clustered. This is to allow the grid
!!    to refine and make it possible to spread out the clustering.
!!    The argument "partPosInitialized" indicates to the caller whether the whole
!!    file was read in. The routine keeps track of the point up to 
!!    which the file was read, and it is called again, it starts 
!!    reading the next values from the current pointer, and not the 
!!    beginning.
!!
!! ARGUMENTS
!!
!!  partPosInitialized : boolean indicating whether particles positions were 
!!            successfully initialized. This is not really relevant
!!            for this version of the routine
!! updateRefine : is true if the routine wished to retain the already
!!                initialized particles instead of reinitializing them
!!                as the grid refine.
!!
!!
!!***


subroutine Particles_initPositions (partPosInitialized,updateRefine)

  use Particles_data, ONLY:  pt_numLocal, particles, pt_maxPerProc, &
       pt_numAtOnce,pt_posInitialized, pt_meshComm
  
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkBoundBox
  use pt_interface, ONLY : pt_numToRead,  pt_initIfInBlock, useParticles

  implicit none
#include "constants.h"

#include "Flash.h"
  include "Flash_mpi.h"

  logical, intent(INOUT) :: partPosInitialized
  logical, intent(OUT) :: updateRefine

  integer       :: i, j, k, n, m, ierr

  integer       :: blockID

  integer       :: blkCount, maxNumParticles
  integer,dimension(MAXBLOCKS) :: blkList
  real, dimension(LOW:HIGH,MDIM) :: boundBox
  logical :: notFound

!----------------------------------------------------------------------

  !       Initialization now done in Particles_init.
  
  !        Particle slot number
  
  if(.not.useParticles) then
     partPosInitialized = .true.
     updateRefine=.false.
  end if
  
  updateRefine=.true.
  if(partPosInitialized) return
  
  call MPI_AllReduce(pt_numLocal,maxNumParticles, 1, MPI_INTEGER, MPI_MAX, pt_meshComm, ierr)  
  numToGet=pt_maxPerProc-maxNumParticles
  
  call Grid_getListOfBlocks(LEAF,blkList,blkCount)
  
  call pt_readParticles(particles,pt_numLocal,numToRead,partPosInitialized)
  j=pt_numLocal
  blkInd=0
  do i=1,numToRead
     j=j+1
     notFound=.true.
     do while(notFound)
        blkInd=blkInd+1
        blockID=blkList(blkInd)
        call Grid_getBlkBoundBox(blockID,boundBox)
        notFound=(boundBox(LOW,IAXIS)>particles(POSX_PART_PROP,j)).or.&
             (boundBox(HIGH,IAXIS)<particles(POSX_PART_PROP,j))
        if(NDIM>1) &
             notFound=(boundBox(LOW,JAXIS)>particles(POSY_PART_PROP,j)).or.&
             (boundBox(HIGH,JAXIS)<particles(POSY_PART_PROP,j)).or.&
             notFound
        
        if(NDIM>2) &:
        notFound=(boundBox(LOW,KAXIS)>particles(POSZ_PART_PROP,j)).or.&
             (boundBox(HIGH,KAXIS)<particles(POSZ_PART_PROP,j)).or.&
             notFound
        if(.not.notFound) then
           pt_numLocal=pt_numLocal+1
           particles(:,pt_numLocal)=particles(:,j)
           particles(BLK_PART_PROP,pt_numLocal)=blockID
        end if
     end do
  end do
  return
  
end subroutine Particles_initPositions


