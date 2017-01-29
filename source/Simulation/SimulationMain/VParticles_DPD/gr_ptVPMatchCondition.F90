!!****if* source/Simulation/SimulationMain/VParticles_DPD/gr_ptVPMatchCondition
!!
!! NAME
!!
!!  gr_ptVPMatchCondition
!!
!! SYNOPSIS
!!
!!  gr_ptVPMatchCondition(real(INOUT)  :: particle(:),
!!                    integer(IN)   :: blockID,
!!                    logical(IN)   :: newBlkID,
!!                   integer(INOUT) :: vpartCount, 
!!                      real(OUT)   :: vparticles)
!!
!! DESCRIPTION
!!     
!!    This routine is used to find if the condition for generating
!!    virtual copies of a particle is met. For this specific implementation
!!    the condition is met if the particle is in the cell adjacent to the
!!    block boundary, whether it is an interior cell, or a guard cell
!!
!! ARGUMENTS
!!
!!     pos -  the position coordinates
!!     blockID   - the block to which the current particle is attached
!!     newBlkID  - indicates if the blockid is new, and therefore
!!                 boundblock and surrblks are needed
!!     vpartCount - number of copies of particles to be generated
!!                  it comes in with maximum value which is 7, and
!!                  returns with the actual number, which can be 0
!!                  if no copies are needed
!!     vparticles - each copy of the particle goes to one neighboring
!!                  block of the current block. The first three indices
!!                  in the array indicate the faces which the neighbor 
!!                  shares with the current block, and the last index
!!                  is set to 1 if that neighbor should get the real
!!                  particle, or 0 if the neighbor gets a virtual copy
!!                  if the real particle is not leaving the block yet,
!!                  but needs copies because it is close to the block
!!                  boundary then all vparticles entries will have a 0
!!                  in the last index.
!! 
!!***
!#define DEBUG_VPARTICLES

#include "Flash.h"
#include "constants.h"

subroutine gr_ptVPMatchCondition(pos, blockID, newBlkID, &
     vpartCount, vparticles)

  use gr_ptVPData, ONLY : gr_ptVPBndBox, gr_ptVPDeltas
  use Simulation_data, ONLY: rc
  implicit none

  real,dimension(MDIM) :: pos
  integer, intent(IN) :: blockID
  logical,intent(IN) :: newBlkID
  integer, intent(INOUT) :: vpartCount
  integer, dimension(4,vpartCount),intent(OUT) :: vparticles

  integer,dimension(MDIM) :: faces,leavingFace
  real :: dist
  integer :: count
  logical :: ax(MDIM),iandj,iandk,jandk,ijk
  integer :: i

  if(newBlkID) then
     call Grid_getBlkBoundBox(blockID,gr_ptVPBndBox)
     call Grid_getDeltas(blockID,gr_ptVPDeltas)
  end if
  count = 0
  vparticles(IAXIS:KAXIS,:)=CENTER
  vparticles(VP_LEAVE,:)=0
  
  faces=CENTER
  leavingFace=CENTER
  
  ax(:) = .false.     !! ax initialized
  !write(*,*)'rc=',rc
 
  do i = 1,NDIM
     ! write(*,*) 'BLKBOUND',gr_ptVPBndBox(:,i)
     dist=pos(i)-gr_ptVPBndBox(LOW,i)
     
     !if(gr_ptVPDeltas(i)>abs(dist))then
     if (rc > abs(dist)) then
        ax(i)=.true.
        faces(i)=LEFT_EDGE
        count=count+1
        vparticles(i,count)=faces(i)
        if(dist<0) then
           leavingFace(i)=LEFT_EDGE
        end if
     end if
     
     dist=gr_ptVPBndBox(HIGH,i)-pos(i)
        !if(gr_ptVPDeltas(i)>abs(dist)) then
     if (rc > abs(dist)) then
        ax(i)=.true.
        faces(i)=RIGHT_EDGE
        count=count+1
        vparticles(i,count)=faces(i)
        if(dist<0) then
           leavingFace(i)=RIGHT_EDGE
        end if
     end if
     !end if
  end do

  iandj=ax(IAXIS).and.ax(JAXIS)
  iandk=ax(IAXIS).and.ax(KAXIS)
  jandk=ax(JAXIS).and.ax(KAXIS)
  ijk=ax(IAXIS).and.ax(JAXIS).and.ax(KAXIS)
  if(iandj)then
     count=count+1
     vparticles(IAXIS,count)=faces(IAXIS)
     vparticles(JAXIS,count)=faces(JAXIS)
  end if
  if(iandk)then
     count=count+1
     vparticles(IAXIS,count)=faces(IAXIS)
     vparticles(KAXIS,count)=faces(KAXIS)
  end if
  if(jandk)then
     count=count+1
     vparticles(JAXIS,count)=faces(JAXIS)
     vparticles(KAXIS,count)=faces(KAXIS)
  end if
  if(ijk)then
     count=count+1
     vparticles(IAXIS:KAXIS,count)=faces(IAXIS:KAXIS)
  end if
  vpartCount=count

  if(count>0) then
     do i = 1,count
        if((vparticles(IAXIS,i)==Leavingface(IAXIS))&
             .and.(vparticles(JAXIS,i)==Leavingface(JAXIS))&
             .and.(vparticles(KAXIS,i)==Leavingface(KAXIS)))&
             vparticles(VP_LEAVE,i)=1
     end do
  end if
end subroutine gr_ptVPMatchCondition
