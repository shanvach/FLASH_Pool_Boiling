!!****if* source/Simulation/SimulationMain/VParticles_DPD/gr_ptVPGenerate
!!
!! NAME
!!
!!  gr_ptVPGenerate
!!
!! SYNOPSIS
!!
!!  gr_ptVPGenerate(real(INOUT)   :: particle(:),
!!                    integer(IN)   :: propCount,
!!                   integer(OUT)   :: destCount, 
!!                      real(OUT)   :: destParticles,
!!                    integer(IN)   :: blockID,
!!                    logical(IN)   :: newBlkID)
!!
!! DESCRIPTION
!!     
!!    This routine is used in moving the non stationaly data elements 
!!    associated with structures like particles and ray, when a data element
!!    moves off a block without re-gridding. Here every element currently 
!!    on the processor is examined to see if it still belongs to the same block.
!!    If it does not, it is further examimned to see if it has moved out of the physical boundary.
!!    If is out of physical boundary, it may either leave the domain, stay on the same block
!!    or be moved  to destBuf, which holds elements to be passed to the next processor, depending
!!    on the boundary conditions. If it is still in the physical domain, it may have
!!    moved to another block on the same processor, in which case only its BLK
!!    needs to change, otherwise it is moved to destBuf.
!!
!! ARGUMENTS
!!
!!     particle -  single element with all its attributes
!!                           values
!!     propCount - number of element attributes
!!
!!     destCount - the count of virtual particles created
!!     destParticles - the newly created virtual particles
!!     blockID   - the block to which the current particle is attached
!!     newBlkID  - indicates if the blockid is new, and therefore
!!                 boundblock and surrblks are needed
!! 
!!***

#ifdef DEBUG_ALL
!#define DEBUG_PARTICLES
#endif
!#define DEBUG_VPARTICLES

subroutine gr_ptVPGenerate(particle,propCount,destCount,destParticles,&
                           blockID,newBlkID)
#include "constants.h"
#include "Flash.h"
#include "gr_ptMapToMesh.h"

#ifdef  DEBUG_VPARTICLES
  use Driver_data,ONLY: dr_globalMe
#endif
  use gr_ptData, ONLY : gr_ptProc,gr_ptBlk,gr_ptPosx,gr_ptPosy,gr_ptPosz
  use gr_ptVPData, ONLY : gr_ptVPMaxCount
  use Grid_data,ONLY : gr_globalDomain, gr_domainBC
  use Simulation_data,ONLY: domainsize
  !-----  Argument list ------
  implicit none
  integer, intent(IN) :: propCount
  real,dimension(propCount),intent(INOUT):: particle
  integer,intent(OUT) :: destCount
  real,dimension(propCount, gr_ptVPMaxCount), intent(OUT) :: destParticles
  integer,intent(IN) :: blockID
  logical, intent(IN) :: newBlkID
  !-------------------------------------
  
  ! local variables
  ! integer :: blockID, procID
  integer,parameter :: MAXCOUNT=14
  integer,dimension(BLKNO:PROCNO,MAXCOUNT) :: neghID
  integer,dimension(BLKNO:PROCNO)::destNeghID
  integer,dimension(BLKID:REFLEVELDIF,ABSMAXNEGH):: negh
  integer,dimension(MDIM,ABSMAXNEGH) :: neghCornerID
  integer :: numNegh
  real,dimension(MDIM) ::  pos
  !real,dimension(MDIM) :: domainsize
  integer,dimension(MDIM) ::  onBoundary
  logical :: leftDomain,boundaryside
  integer :: vpartCount,side
  integer,dimension(IAXIS:VP_LEAVE,MAXCOUNT)::vparticles
  integer i,j
  ! ------ end of local variables -----------------

  ! Initialization of the number of virtual particles that will be created
  destCount = 0 
  !domainsize(IAXIS:KAXIS)=0.;


  ! Check if a real particle has left the domain
  call gr_ptVPBC(particle,propCount, leftDomain, onBoundary)

#ifdef DEBUG_VPARTICLES
  !write(*,*)'pos=',particle(gr_ptPosx:gr_ptPosz),'left=',leftDomain,'onBoundary',onBoundary
  if (sum(onBoundary)>-1368) then
     do i=1,NDIM
        side=onBoundary(i);
        if ((side==1)) then
           if (gr_domainBC(side,i)==PERIODIC) &
                print *,'This particle is on a left periodic boundary'
           if (gr_domainBC(side,i)==REFLECTING) &
                print *,'This particle is on a left reflecting boundary'
        elseif (side==2) then 
           if (gr_domainBC(side,i)==PERIODIC) &
                print *,'This particle is on a right periodic boundary'
           if (gr_domainBC(side,i)==REFLECTING) &
                print *,'This particle is on a right reflecting boundary'
        end if
     end do     
  end if
#endif

  if(leftDomain)then !! nothing needs to be done except 
                     !!letting the particle cease to exist 
     !this section needs to be updated to change the positions of the particle in case of PERIODIC BC.
     write(*,*) 'A particle left the domain'
     stop
     particle(gr_ptBlk)=LOST
     destCount=-1
  else
      
     vpartCount=MAXCOUNT
     pos(1:MDIM)=particle(gr_ptPosx:gr_ptPosz)

     !write(*,*)'proc=',int(particle(gr_ptProc)),'vpart before=',vpartCount
     call gr_ptVPMatchCondition(pos,blockID,newBlkID,vpartCount,vparticles)
     !write(*,*)'proc=',int(particle(gr_ptProc)),'vpart',vpartCount,'needs handling'
     if(vpartCount>0) then !! the particle needs handling 
        !write(*,*)'vp=',vparticles(VP_LEAVE,1:vpartCount)
        do i = 1,vpartCount
            
           if(vparticles(VP_LEAVE,i)==1) then
              !! This part is necessary when the particle is in the center of a face
              !! and the neighbor along that face is refined. Then potentially there
              !! there is a possibility of four copies of the particles, and if it is
              !! leaving the current block we need to determine which of the four 
              !! neghbors is the destination
             
              call gr_findNeghID(blockID,pos,vparticles(IAXIS:KAXIS,i),destNeghID)
              !write(*,*)'particle leaving'
              !write(*,*) 'Destination BLK and PROC',destNeghID
           end if

           !write(*,*) 'Destination BLK and PROC',destNeghID
           !! And now we find out if we missed any neghbors is call to match condition
           !! because of fine-coarse boundary at this neighbor
           !write(*,*) 'Checking for particle pos',  pos(1:NDIM),'TAG',particle(TAG_PART_PROP)
 !          call gr_ptFindNegh(blockID,vparticles(IAXIS:KAXIS,i),neghid,neghCornerID,numNegh)



#ifdef DEBUG_VPARTICLES
           if (dr_globalMe==2)then
              write(*,*)'vp I:K=',vparticles(IAXIS:KAXIS,i),'proc',dr_globalMe
              write(*,*)'numNegh',numNegh,'neghid',neghid(:,numNegh)
              write(*,*)'pos',pos(1:3)
           end if
#endif

           
           !! If there was a fine-coarse boundary and the particle is sitting close to
           !! to the middle of that boundary numNegh will be greater than 1
           !! destCount = 0      
           do j = 1,numNegh
              destCount=destCount+1
              destParticles(:,destCount)=particle(:)
              if     (onBoundary(IAXIS)==vparticles(IAXIS,i)) then
                 if ((onBoundary(IAXIS)==LEFT_EDGE).and.&
                      (gr_domainBC(LOW,IAXIS)==PERIODIC)) then
                 destParticles(gr_ptposx,destCount)=particle(gr_ptPosx)+domainsize(IAXIS);
                 elseif ((onBoundary(IAXIS)==RIGHT_EDGE).and.&
                      (gr_domainBC(HIGH,IAXIS)==PERIODIC)) then
                    destParticles(gr_ptposx,destCount)=particle(gr_ptPosx)-domainsize(IAXIS);
                 end if

#ifdef DEBUG_PERIODIC
                 write(*,*)'Old x pos',particle(gr_ptPosx), &
                      'New x pos',destParticles(gr_ptposx,destCount)
#endif
              endif
              if (onBoundary(JAXIS)==vparticles(JAXIS,i)) then
                  if ((onBoundary(JAXIS)==LEFT_EDGE).and.&
                      (gr_domainBC(LOW,JAXIS)==PERIODIC)) then
                 destParticles(gr_ptposy,destCount)=particle(gr_ptPosy)+ domainsize(JAXIS);
                 elseif  ((onBoundary(JAXIS)==RIGHT_EDGE).and.&
                      (gr_domainBC(HIGH,JAXIS)==PERIODIC)) then
                     destParticles(gr_ptposy,destCount)=particle(gr_ptPosy)- domainsize(JAXIS);
                 end if
#ifdef DEBUG_PERIODIC
                 write(*,*)'Old y pos',particle(gr_ptPosy), &
                'New y pos',destParticles(gr_ptposy,destCount)
#endif
              endif
#if NDIM==3
              if (onBoundary(KAXIS)==vparticles(KAXIS,i)) then
                  if ((onBoundary(KAXIS)==LEFT_EDGE).and.&
                      (gr_domainBC(LOW,KAXIS)==PERIODIC)) then
                 destParticles(gr_ptposz,destCount)=particle(gr_ptPosz)+ domainsize(KAXIS);
                 elseif ((onBoundary(KAXIS)==RIGHT_EDGE).and.&
                      (gr_domainBC(HIGH,KAXIS)==PERIODIC)) then
                     destParticles(gr_ptposz,destCount)=particle(gr_ptPosz)- domainsize(KAXIS);
                  end if
#ifdef DEBUG_PERIODIC
                  write(*,*)'Old z pos',particle(gr_ptPosz), &
                       'New z pos',destParticles(gr_ptposz,destCount)
#endif
              end if
#endif
              
              destParticles(PROC_PART_PROP,destCount)=neghid(PROCNO,j)
              destParticles(BLK_PART_PROP,destCount)=neghid(BLKNO,j)
              !! If the neghid found by gr_ptFindNegh matches with
              !! that found by gr_findNeghID then this is the real
              !! copy of the particle
              if((neghid(BLKNO,j)==destNeghID(BLKNO)).and.((neghid(PROCNO,j)==destNeghID(PROCNO))))then
                 particle(TAG_PART_PROP)=-particle(TAG_PART_PROP)
                 
#ifdef DEBUG_VPARTICLES
                 write(*,*) 'Real particle i=',i,'new TAG=',int(particle(TAG_PART_PROP)),'on proc',int(particle(PROC_PART_PROP))
#endif
                 
              else
                 destParticles(TAG_PART_PROP,destCount)=-particle(TAG_PART_PROP)
                 
#ifdef DEBUG_VPARTICLES
                 write(*,*) destCount,'on proc=',int(destParticles(PROC_PART_PROP,destCount)),' TAG=', & 
                      destParticles(TAG_PART_PROP,destCount)
#endif
              end if
           end do
           
        end do
      
     end if
  end if
  !write(*,*)'Final no. of virtual particles created is= ',destcount
end subroutine gr_ptVPGenerate
