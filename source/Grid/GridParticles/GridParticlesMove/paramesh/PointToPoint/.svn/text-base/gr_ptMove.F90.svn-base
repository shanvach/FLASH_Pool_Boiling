!!****if* source/Grid/GridParticles/GridParticlesMove/paramesh/PointToPoint/gr_ptMove
!!
!! NAME
!!  gr_ptMove
!!
!! SYNOPSIS
!!
!!  gr_ptMove(real(INOUT)    :: databuf(:,:),
!!            integer(IN)    :: propCount 
!!            integer(INOUT) :: localCount,
!!            integer(IN)    :: maxCount,
!!            logical(INOUT) :: moveDone)
!!                    
!!  
!! DESCRIPTION 
!!  
!!  This routine moves the particles data to their destination blocks
!!  after time integration
!!
!!  With time integration only a small fraction of particles move out
!!  of a block at any timestep. However, all particles must be examined
!!  to determine if they moved out of their curret block. With refinement
!!  all particles of one block move together to a new block. The logistics
!!  of moving the data between processors is the same in both situations.
!!  Therefore this routine can be used in both modes. 
!! 
!!  Overview of algorithm
!!
!!  * Call a local matching routine to find which of the displaced
!!    particles have moved to another block within the same processor
!!    and move them there.  For local move
!!    from one block to another, the particle just needs the block number
!!    field changed in its data structure, there is no actual data movement.
!!    Note that local matching algorithms are different for the two modes.
!!
!!  * Sort the particles in the send buffer by their destination processor
!!
!!  * Do a global communication to let each processor know the number of other
!!    processors that will be sending particles to it.
!!
!!  * Post non-blocking receives equal to the number of processors found in the 
!!    previous step, and send particles in the send buffer to their destination
!!  
!!  * Place the received particles in their appropriate blocks
!!
!!
!! ARGUMENTS 
!!
!!  databuf : List of particles. It is a 2D real array, the first dimension
!!              represents particle's properties, and second dimension is 
!!              index to particles.
!!
!! propCount : number of properties for this datastructure 
!!
!!  localCount : While coming in it contains the current 
!!               number of elements in the data structure mapped to
!!               this processor. After all the data structure 
!!               movement, the number might change, 
!!               and the new value is put back into it
!!  maxCount : This is parameter determined at runtime, 
!!             and is the maximum count of elements 
!!             that a simulation expects to have. 
!!             All the arrays  are allocated based on this number
!!  moveDone  : a logical argument to indicate whether the particles
!!              have been moved to their destination
!!
!!
!! NOTES
!!   
!!
!! SEE ALSO
!!
!!
!!
!!
!!***

#define GLOBAL

subroutine gr_ptMove(databuf,propCount, localCount,globalCount, moveDone)

  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_data, ONLY : gr_useParticles, gr_meshNumProcs, gr_meshMe, gr_meshComm, &
       gr_useEnergyDeposition
  use gr_ptData, ONLY : gr_ptBlkList,gr_ptBlkCount,&
       gr_ptDestBuf, gr_ptSourceBuf, gr_ptBlk, gr_ptProc
  use ut_sortInterface, ONLY : ut_sortOnProcs
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer,intent(INOUT) :: localCount
  integer,intent(IN) :: globalCount, propCount

  real, dimension(propCount, globalCount),intent(INOUT) :: databuf
  logical, intent(INOUT) :: moveDone

  integer :: numDest,blockID,ierr

  integer, dimension(gr_meshNumProcs) :: perProc,toProcs, fromProcs,&
       fromCount, maxCount
  integer, allocatable, dimension(:,:) :: status
  integer, allocatable, dimension(:) :: req
  integer :: sendCount,recvCount, i,j,k
  integer :: countAndProc,bufSize
  integer :: gettingFrom, sendingTo
  integer :: logUnit
  logical, parameter :: logUnitLocal=.true.
  integer, parameter :: tag_for_count=23
  integer, parameter :: tag_for_data=33

  !if we have turned off particles then return
  
  if(.not. (gr_useParticles .or. gr_useEnergyDeposition)) return

  call Timers_start("gr_ptMove")

  numDest=0

  
  !! make sure that the space in the particles array that is not populated
  !! by existing particles doesn't have valid block number 
  databuf(gr_ptBlk,localCount+1:globalCount)=NONEXISTENT

  !! find particles that have moved off block. The moved particles are
  !! returned in gr_ptDest Buf
  call gr_ensureValidNeighborInfo(0)
  call gr_ptMoveOffBlk(databuf,propCount,localCount,globalCount,gr_ptDestBuf,numDest)

  if(gr_meshNumProcs==1) then
     call Timers_stop("gr_ptMove")
     moveDone = .true.
     ! Updating the dataBuf and localCount with the virtual copies created and stored in gr_DestBuf
     dataBuf(:,localCount+1:localCount+numDest)=gr_ptDestBuf(:,1:numDest)
     
     localCount = localCount+numDest
    
     return
  end if
  
  
  
  fromCount=0
  !! Sort the particles to be moved based upon the processor number of 
  !! their destination
  call ut_sortOnProcs(numDest,propCount, gr_ptProc, gr_meshNumProcs,&
        gr_ptDestBuf,gr_ptSourceBuf, &
       perProc, toProcs,sendingTo)

  !! share that information with all other processors 
  call MPI_Allreduce(toProcs,fromProcs,gr_meshNumProcs,FLASH_INTEGER,MPI_SUM,gr_meshComm,ierr)
  !! At this point all the processors know how many processors they should
  !! expect to get the data from. That information is stored in fromProcs
  !! the +1 is there because proc numbers start from zero.
  
  !! The communication itself happens in two steps. In the first each
  !! processor gets the number of particles to expect from all the processor
  !! that are sending to it, and in the next step the actual transmission of
  !! particles takes place



  gettingFrom=fromProcs(gr_meshMe+1)

  allocate(status(MPI_STATUS_SIZE,gettingFrom))
  allocate(req(gettingFrom))
  req(:)=0
  status=0
  recvCount = 1

#ifdef GLOBAL
  call MPI_ALLReduce(perProc,maxCount, gr_meshNumProcs, FLASH_INTEGER,MPI_MAX,gr_meshComm,ierr)
 recvCount=maxCount(gr_meshMe+1)
#else
  !! Post all non-blocking receives necessary
  if(gettingFrom>0) then
     do i = 1,gettingFrom
        call MPI_Irecv(fromProcs(i),recvCount,FLASH_INTEGER,MPI_ANY_SOURCE,&
             tag_for_count,gr_meshComm,req(i),ierr)
     end do
  end if

  !! schedule all sends in blocking mode
  if(sendingTo>0) then
     j=0
     do i=1,sendingTo
        do while(perProc(j+1)==0)
           j=j+1
        end do
        
        countAndProc=perProc(j+1)
        sendCount=1
        call MPI_Send(countAndProc,sendCount,FLASH_INTEGER,j,tag_for_count,gr_meshComm,ierr)
        j=j+1
        
     end do

  end if

  !! This is where we wait for the all the posted receives to be done,
  !! find the maximum number of particles mype is expected to receive
  !! and allocate a buffer large enough to receive it.

  if(gettingFrom>0) then
     call MPI_Waitall(gettingFrom,req,status,ierr)
     recvCount=(maxval(fromProcs(1:gettingFrom)))
  end if
#endif
  if(gettingFrom>0) then
     if ((gettingFrom*recvCount) > ubound(gr_ptSourceBuf,2)) then
        print *, "Overflow! Mesh PE", gr_meshMe, &
             ", gettingFrom", gettingFrom, &
             ", recvCount", recvCount, &
             ", max particles", ubound(gr_ptSourceBuf,2)
        call Driver_abortFlash("[gr_ptMove]: Insufficient space "//&
             "in particles communication buffer: "//&
             "increase pt_maxPerProc")
     end if
     gr_ptSourceBuf(:,:)=NONEXISTENT
     bufSize = recvCount*propCount
     
     req(:)=0
     status=0
     j=1

     !! Here the receives are posted to get the particles
     do i =1,gettingFrom
        call MPI_Irecv(gr_ptSourceBuf(1,j),bufSize,&
             FLASH_REAL,MPI_ANY_SOURCE,tag_for_data,gr_meshComm,&
             req(i),ierr)
        j=j+recvCount
     end do
     
  end if

  if(sendingTo>0) then
     j=0
     k=1
     
     do i=1,sendingTo
        do while(perProc(j+1)==0)
           j=j+1
        end do
        sendCount=perProc(j+1)
        bufSize=sendCount*propCount
        call MPI_Send(gr_ptDestBuf(1,k),bufSize,FLASH_REAL,j,tag_for_data,gr_meshComm,ierr)
        j=j+1
        k=k+sendCount
        
     end do
  end if
  

  if(gettingFrom>0) then
     call MPI_Waitall(gettingFrom,req,status,ierr)
     do i = 1,recvCount*gettingFrom
        
        blockID=int(gr_ptSourceBuf(gr_ptBlk,i))
        if(blockID/=NONEXISTENT) then
           localCount=localCount+1
           if (localCount > globalCount) then
              call Driver_abortFlash("[gr_ptMove]: Insufficient space "//&
                   "in particles array: increase pt_maxPerProc")
           end if
           databuf(:,localCount)=gr_ptSourceBuf(:,i)
        end if
     end do
     deallocate(req)
     deallocate(status)
  end if

  call Timers_stop("gr_ptMove")
  moveDone = .true.
end subroutine gr_ptMove


