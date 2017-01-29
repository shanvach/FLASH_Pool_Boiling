!!****if* source/Grid/GridParticles/GridParticlesMove/UG/Directional/gr_ptMoveOffProc
!!
!! NAME
!!
!!  gr_ptMoveOffProc
!!
!! SYNOPSIS
!!
!!   gr_ptMoveOffProc(integer(IN)   :: face, 
!!                   integer(IN)   :: axis,
!!                   integer(IN)   :: index, 
!!                   integer(IN)   :: propCount,
!!                   integer(IN)   :: maxPerProc,
!!                   logical(IN)   :: boundary,
!!                   integer(IN)   :: lnegh,
!!                   integer(IN)   :: rnegh,
!!                   real(IN)      :: corner(LOW:HIGH),
!!                   integer(inout)  :: localNum,
!!                   real,(inout)  :: particles(:,:))
!!
!! DESCRIPTION
!!     
!!    This routine is used in moving the particles data when a particle moves
!!    off a block because of time integration. The routine examines a single face
!!    of the block at a time. Here every particle currently 
!!    on the processor is examined to see if it has moved out of the block across
!!    the specific face and axis under consideration.
!!    If the face was along physical boundary, then it may either leave the domain, 
!!    stay on the same block or be moved  to gr_ptDestBuf, 
!!    which holds particles to be passed to the next processor, depending
!!    on the boundary conditions. Once all particles have been examined, the contents
!!    of gr_ptDestBuf are sent to the neighbor along the face, and the contents
!!    of opposite face neighbors gr_ptDestBuf are received. The newly received particles
!!    are added to the particles data structure. Note that it is assumed that no particle
!!    can move across a full block in a single time step.
!!
!! ARGUMENTS
!!
!!     face - valid values are LOW and HIGH, indicating whether it is lower face or
!!            the upperface along the specified dimension
!!     axis - the physical dimension under consideration
!!     index - index into the particles data structure indicating the particle position 
!!             specified axis
!!     propCount - the count of particle properties
!!     maxPerProc - the maximum number of particles allowed on one processor
!!     boundary - indicates whether the face under consideration is on physical boundary
!!     lnegh     - neghbor to send data to
!!     rnegh     - neghbor to receive data from
!!     corner   - the coordinate of the lower upper face along the axis
!!     localNum -   The number of particles in the particle array that
!!                           are stored on this processor.
!!     particles -           A 2 dimensional array of particles and their property
!!                           values
!!
!!***

subroutine gr_ptMoveOffProc(face,axis,index, propCount, maxPerProc,boundary,&
     lnegh,rnegh,corner,localNum,particles, gr_ptVPDistFromFace)
#include "constants.h"
#include "Flash.h"

  use Grid_data, ONLY : gr_meshMe,gr_axisComm,gr_domainBC
  use gr_ptData, ONLY : gr_ptDestBuf, gr_ptSourceBuf,gr_ptBlk,gr_pttag
  use gr_ptInterface, ONLY : gr_ptOneFaceBC
  use Driver_interface, ONLY: Driver_abortFlash
  implicit none
  include "Flash_mpi.h"


  integer, intent(IN) :: face, axis, index, propCount, maxPerProc
  logical,intent(IN) :: boundary
  integer, intent(IN) :: lnegh,rnegh
  real,dimension(LOW:HIGH),intent(IN) :: corner
  integer,intent(INOUT) :: localNum
  real,dimension(propCount,maxPerProc),intent(INOUT) :: particles
  real, optional, intent(IN) :: gr_ptVPDistFromFace
!!  real, allocatable, dimension(:,:) :: sendBuf,recvBuf

  integer :: j, k, numDest, lostParticles, count,VPcount

  integer :: pid, pend,bufsize, pos
  integer :: recvCount, tag, ierr
  integer,dimension(MPI_STATUS_SIZE) :: status
  real :: nothing, dist_from_face, newPos, ptTag, dist(MDIM)
  logical :: moved, leaving, inVPZone,VP, VP_ptag, hiFace 
  logical :: hiF_VPtag, lcyc
  integer :: blkID

  numDest = 0
  lostParticles=0
  count=0
  tag=20*axis
  nothing = 0.0

  if(face>HIGH)call Driver_abortFlash("Grid_moveParticles: face value is not LOW/HIGH")
  !! find out if working on the high side of the face
  hiFace   = (face .eq. HIGH)
  pend=localNum
  j=1
  VPcount=pend
  do pid = 1,localNum
     lostParticles=0
     VP_ptag = (int(particles(gr_ptTag,j)).lt.0)
     hiF_VPtag = hiFace .and. VP_ptag
      
     blkID=int(particles(gr_ptBlk,j))
     if(blkID.le.0) cycle  !! verify that it is a valid
                           !! particle before processing it
                           !! compute distances from the face 

     dist_from_face = particles(index,j)-corner(face)

!!$     lcyc = .false.
!!$     if (hiF_VPtag) then
!!$        ! For the x axis all Virtual particles on hiFace are cycled.
!!$        select case(axis)
!!$        case(IAXIS) 
!!$           lcyc =.true.
!!$        ! Other axis cases are cycled if the Virtual particle lies out of bounds.
!!$        case default
!!$           if (dist_from_face .gt. 0.) lcyc =.true. 
!!$        end select
!!$     endif
!!$     if (lcyc) cycle

     if (hiF_VPtag .and. (dist_from_face .gt. 0.)) cycle


     !! check if particles could be generating VPS, and if
     !! so whether it is staying in the block or leaving it
     inVPZone=abs(dist_from_face).le.gr_ptVPDistFromFace         
           
     if(inVPZone) then 
        !! determine whether the particle is leaving the current block
        leaving = (inVPZone.and.(((face==LOW ).and.(dist_from_face<0.)).or.&
                                     ((face==HIGH).and.(dist_from_face>0.))))
        !! save tag for use in determining if it needs to be negated later
        ptTag=particles(gr_ptTag,j)
              
        !! if the block is on physical boundary then periodic boundary conditions
        !! need to be managed whether the particles is leaving the block or not,
        !! for every other type of boundary conditions the boundary conditions have
        !! to be managed only if the particles is leaving the block
        VP=.true.
        if(boundary) then
           VP=.false.
           call gr_ptOneFaceBC(particles(:,j),propCount, axis,face,blkID, lostParticles,&
                                  leaving, VP, pos, newPos)
           !! At this point, either the particle will have left the domain
           !! in which case there is nothing further to be done, or if reflective
           !! BC then it should be back in the domain and again does not need
           !! anything further done. Only if it has changed pos that lies outside
           !! the box do we need to create virtual particles, VP will set to true
           !! and newPos will have the new position
        end if
        if(VP) then
           numDest=numDest+1
           gr_ptDestBuf(:,numDest)=particles(:,j)
           !! if the particle is going across a corner then it may have already
           !! generated one virtual particle along an earlier dimension, here
           !! we are generating a virtual particle off of a virtual particle
           !! the tag here will already be negative, so doesn't need any 
           !! special processing
           !! If the particle was real then we have to figure out
           !! which one is the real particle and which one is virtual
           if(leaving) then
              if(ptTag>0) particles(gr_ptTag,j)=-ptTag !! if real pos of the particle is 
                       !! is outside the block then in the block the particle is virtual
           else
              if(ptTag>0) gr_ptDestBuf(gr_ptTag,numDest)=-ptTag 
                       !! otherwise the copy going to the
                       !! neghboring block is virtual
           end if
           !! if the crossing happened on the boundary then the copy going to the
           !! neighbor block also needs the new position
           if(boundary)gr_ptDestBuf(pos,numDest)=newPos
        end if
     end if
     if(lostParticles>0) then
        particles(:,j)=particles(:,pend)
        pend=pend-1
     else
        j=j+1
     end if
  end do

  localNum=pend
  count=numDest*propCount
  if(count==0)then   !! if there were no particles to send to the neighbor, 
     count=1         !! send a single data item so sendreceive doesn't hang
     gr_ptDestBuf(1,1)=nothing
  end if

  bufsize=propCount*maxPerProc

  recvCount=bufSize


  call MPI_Sendrecv(gr_ptDestBuf,count,FLASH_REAL,lnegh,tag,&  !! send to left neghbor, receive from 
       gr_ptSourceBuf,recvCount,FLASH_REAL,rnegh,tag,&         !! the right one
       gr_axisComm(axis),status,ierr)
  call MPI_Get_count(status,FLASH_REAL,recvCount,ierr)         !! find out how many particle received

  if(recvCount>1)then       !! if some particles were received, then append them to
     k=recvCount/propCount
     do j = 1,k
        gr_ptSourceBuf(PROC_PART_PROP,j) = gr_meshMe
        particles(:,localNum+j)=gr_ptSourceBuf(:,j)
     end do
     !write(*,*) 'REAL - VIRTUAL PARTS=',gr_meshMe,axis,localnum,k
     localNum=localNum+k
  end if
end subroutine gr_ptMoveOffProc
