!!****if* source/Particles/ParticlesInitialization/FromFile/pt_initIfInBlock
!!
!! NAME
!!    pt_initIfInBlock
!!
!! SYNOPSIS
!!
!!    call pt_initIfInBlock(real(IN)       :: boundBox(LOW:HIGH,MDIM),
!!                          integer(INOUT) :: numParticles,
!!                          integer(IN)    :: putStart,
!!                          integer(INOUT) :: getEnd)
!!
!!
!! DESCRIPTION
!!    
!!    Given a set of particles, this routine determines if the particles are
!!    contained within the bounding box provided by the caller. The contained 
!!    particles are moved to the front end of the particles data structure, 
!!    whereas the unfound ones stay at the back end.
!!    
!!
!! ARGUMENTS
!!
!!  boundBox - the coordinates of the local lefthand and upper righthand corner
!!             of a box to be considered by this routine
!!  numParticles - at entry it indicates the number of particles to be processed
!!                 at exit in contains the number of particle contained within
!!                 the specified boundbox
!!  putStart     - the starting index into the particles data structure where space
!!                 is available to start putting in the found particles
!!  getEnd       - getEnd it the index at the end of remaining particles to be processed
!!                 Upon entry it points to the end of the particles data structure
!!                  As particles get found, pointer also moves, always pointing to the 
!!                 last valid unprocessed particle in the data structure.
!!
!!
!!***


subroutine pt_initIfInBlock(boundBox,numParticles,putStart,getEnd)

  use Particles_data, ONLY:  particles
  
  implicit none
#include "constants.h"

#include "Flash.h"

  real, dimension(LOW:HIGH,MDIM), intent(IN) :: boundBox
  integer, intent(IN) :: putStart
  integer, intent(INOUT) :: getEnd,numParticles

  integer       :: getStart
  integer       :: i, j, k, n, m

  logical       :: IsInBlock


  getStart=getEnd-numParticles+1

  j=putStart
  m=getEnd
  n=getStart
  numParticles=0
  do i = getStart,getEnd
     IsInBlock=Particles(POSX_PART_PROP,n)>=boundBox(LOW,IAXIS)
     if(IsInBlock)IsInBlock=Particles(POSX_PART_PROP,n)<=boundBox(HIGH,IAXIS)
     if(NDIM>1) then
        if(IsInBlock)IsInBlock=Particles(POSY_PART_PROP,n)>=boundBox(LOW,JAXIS)
        if(IsInBlock)IsInBlock=Particles(POSY_PART_PROP,n)<=boundBox(HIGH,JAXIS)
     end if
     if(NDIM>2) then
        if(IsInBlock)IsInBlock=Particles(POSZ_PART_PROP,n)>=boundBox(LOW,KAXIS)
        if(IsInBlock)IsInBlock=Particles(POSZ_PART_PROP,n)<=boundBox(HIGH,KAXIS)
     end if
     if(IsInBlock) then
        j=j+1
        Particles(:,j)=Particles(:,n)
!        Particles(BLK_PART_PROP,j)=blockID
        Particles(BLK_PART_PROP,j)=1
        if(m>n)Particles(:,n)=Particles(:,m)
        m=m-1
        numParticles=numParticles+1
     else
        n=n+1
     end if
  end do
  getEnd=m
  
  return
  
  !----------------------------------------------------------------------
  
end subroutine Pt_initIfInBlock


