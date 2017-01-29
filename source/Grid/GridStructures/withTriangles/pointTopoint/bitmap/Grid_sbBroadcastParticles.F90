!!****if* source/Grid/GridStructures/withTriangles/pointTopoint/bitmap/Grid_sbBroadcastParticles
!!
!! NAME
!!  Grid_sbBroadcastParticles
!!
!! SYNOPSIS
!!
!!  Grid_sbBroadcastParticles()
!!  
!! DESCRIPTION 
!!  
!!  The particle broadcast routine
!!
!!  Overview of the algoritm
!!
!!  * Called from Grid_initDomain
!!
!!  * Gets the coordinates of all particles. Calculate the number of particles 
!!  that should be sent from the master to all processors
!!
!!  * Scatter the position of the particles that do not belong to the Master PE
!!    to all processors
!!
!!  * Find which processor the particle belongs to
!!
!!  * Update forces
!!
!!  * Each processor sends updated forces back to master
!!
!! ARGUMENTS 
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine Grid_sbBroadcastParticles()
  use Timers_interface,    ONLY : Timers_start,                  &
       Timers_stop,                   &
       Timers_getSummary
 
!  use gr_sbData, ONLY : gr_sbParticleCount
  
  implicit none
  call Timers_start("body_broadcast_particles")
  call gr_sbGetProcBlock()
  call gr_sbStoreParticlesPerProc()
  call gr_sbSendParticleCount()
  call gr_sbSendParticles()

!  At the start of the computations we dont do IB forcing - MV
!  call gr_sbUpdateForces()
!  call gr_sbSendForces()
!  if (allocated(gr_sbParticleCount)) deallocate(gr_sbParticleCount)

!   use Grid_data, ONLY : gr_meshMe, gr_meshComm, gr_meshNumProcs
!   use Grid_interface, ONLY : Grid_xyzToBlock
!   use bittree, only : amr_identify_block
!   use RuntimeParameters_interface, ONLY : RuntimeParameters_get
!   use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, totalPart, sumPart
!   use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkBoundBox, &
!        Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr
!   use Timers_interface, ONLY : Timers_start, Timers_stop
!   implicit none
!   include "Flash_mpi.h"

!   integer :: i, j, k, blkCount, count, recv, recvBufCount, ierr, blkID, lev, &
!        proc, blk, b, localParticleCount, bufSize, send, recvSize, sendBufCount, max
!   real, dimension(2,MDIM) :: boundBox
!   real, dimension(MDIM) :: particleposn, particleposnSB
!   integer,dimension(MAXBLOCKS) :: listOfBlocks
!   real, dimension(NPART_PROPS,7600) :: SourceBuf
!   integer, dimension(MDIM) :: ijk
!   integer, dimension(0:gr_meshNumProcs-1) :: fromProc
!   integer, allocatable, dimension(:) :: sendreq, req
!   integer, allocatable, dimension(:,:) :: sstatus, rstatus

!   call Timers_start("body_broadcast_particles")

!   allocate(sstatus(MPI_STATUS_SIZE, gr_meshNumProcs))
!   allocate(rstatus(MPI_STATUS_SIZE, gr_meshNumProcs))
!   allocate(req(gr_meshNumProcs))
!   req = 0
!   sstatus = 0
!   recvBufCount = NPART_PROPS
!   sendBufCount = NPART_PROPS
!   send = 0
!   recv = 0
!   recvSize = 0
!   allocate(sendreq(gr_meshNumProcs))
!   sendreq = MPI_REQUEST_NULL

!   call RuntimeParameters_get("lrefine_max", lev)

! !  do b = 1, gr_sbNumBodies
! !     if (gr_sbBodyInfo(b) % myPE == gr_sbBodyInfo(b) % bodyMaster) then
! !        if (b == 1) then
! !           max = gr_sbBodyInfo(1) % totalPart
! !        end if
! !        if (gr_sbBodyInfo(b) % totalPart > max) then
! !           max = gr_sbBodyInfo(b) % totalPart
! !        end if
! !     end if
! !  end do

! !  CALL MPI_ALLREDUCE(MPI_IN_PLACE, max,1,&
! !          FLASH_INTEGER, MPI_SUM,gr_meshComm,ierr)
! !  allocate(SourceBuf(NPART_PROPS,max))

!   do b = 1, gr_sbNumBodies
!      fromProc = 0
!      bufSize = 0
!      SourceBuf = 0
!      send = 0
!         ! If master PE, the Master PE gets the positions and the processors
!         ! of all particles
!      if (gr_sbBodyInfo(b) % myPE == gr_sbBodyInfo(b) % bodyMaster) then
!         localParticleCount = 0
!         j = 0
!            ! get Blocks in Master PE
!         call Grid_getListOfBlocks(LEAF, listOfBlocks, count)

!         totalPart = gr_sbBodyInfo(b) % totalPart

!         do i = 1, totalPart
!            particleposn(IAXIS) = gr_sbBodyInfo(b) % particles(POSX_PART_PROP,i)
!            if(NDIM > 1) then
!               particleposn(JAXIS) = gr_sbBodyInfo(b) % particles(POSY_PART_PROP,i)
!            endif
!            if(NDIM > 2) then
!               particleposn(KAXIS) = gr_sbBodyInfo(b) % particles(POSZ_PART_PROP,i)
!            endif
!            call Grid_xyzToBlock(lev, particleposn, ijk)
!            call amr_identify_block(gr_meshNumProcs, lev, ijk, proc, blk)
!            fromProc(proc) = fromProc(proc) + 1
!            do blkCount = 1, count
!               blkID = listofBlocks(blkCount)
!               call Grid_getBlkBoundBox(blkID, boundBox)
!               if (all(particleposn(1:NDIM) >= boundBox(LOW,1:NDIM) .and. &
!                    particleposn(1:NDIM) < boundBox(HIGH,1:NDIM))) then
! !                 localParticleCount = localParticleCount + 1
!                  gr_sbBodyInfo(b) % particles(BLK_PART_PROP,i) = blkID 
!                  gr_sbBodyInfo(b) % particles(PROC_PART_PROP,i) = gr_sbBodyInfo(b) % myPE
!                  !call Grid_updateSolidBodyForces(blkID, b, i, particleposn)
!                  !particle in Master PE
!                  exit
!               endif
!            enddo
!            ! If not in master PE, no. of particles to be sent by Master PE
!            if (blkCount > count) then
!               bufSize = bufSize + NPART_PROPS 
!               j = j + 1
!               SourceBuf(1:NPART_PROPS,j) = &
!                    gr_sbBodyInfo(b) % particles(1:NPART_PROPS,i)
!            endif
!         enddo
!         if (bufSize > 0) then
!            k = 1
!            recv = 0
!            do i = 0, gr_meshNumProcs - 1
!               if (i /= gr_sbBodyInfo(b) % bodyMaster .and. fromProc(i) > 0) then
!                  recvBufCount = fromProc(i)*NPART_PROPS
!                  recv = recv + 1
!                  call MPI_IRECV(gr_sbBodyInfo(b) % particles(1,k), &
!                       recvBufCount, FLASH_REAL, i, &
!                       b, gr_meshComm, req(recv), ierr)
!                  k = k + fromProc(i)
!               endif
!            enddo
!            !print*,"Proc:",gr_meshMe,"K Value:",k
!         endif
!      endif
     
!      ! Broadcast total particle count (numParticles * totalProps)
!      call MPI_BCAST(bufSize, 1, FLASH_INTEGER, &
!           gr_sbBodyInfo(b) % bodyMaster, MPI_COMM_WORLD, ierr)
!      if (bufSize > 0) then
!         ! Broadcast particle positions to all processors
!         call MPI_BCAST(SourceBuf, bufSize, FLASH_REAL, &
!              gr_sbBodyInfo(b) % bodyMaster, MPI_COMM_WORLD, ierr)
!            ! Broadcast particle count (no of particles)
!         call MPI_BCAST(j, 1, FLASH_INTEGER, gr_sbBodyInfo(b) % bodyMaster, &
!              MPI_COMM_WORLD, ierr)
!      endif

!      if (gr_sbBodyInfo(b) % myPE /= gr_sbBodyInfo(b) % bodyMaster) then  !Other processors
!         localParticleCount = 0
!         if (bufSize > 0) then
!            call Grid_getListOfBlocks(LEAF, listOfBlocks, count) !get blocks  
!            do k = 1, j
!               do blkCount = 1, count
!                  blkID = listOfBlocks(blkCount)                    
!                  call Grid_getBlkBoundBox(blkID, boundBox)
!                  particleposnSB(IAXIS) = SourceBuf(POSX_PART_PROP,k)
!                  if(NDIM > 1) then
!                     particleposnSB(JAXIS) = SourceBuf(POSY_PART_PROP,k)
!                  endif
!                  if(NDIM > 2) then
!                     particleposnSB(KAXIS) = SourceBuf(POSZ_PART_PROP,k)
!                  endif
                 
!                  if (all(particleposnSB(1:NDIM) >= boundBox(LOW,1:NDIM) .and. &
!                       particleposnSB(1:NDIM) < boundBox(HIGH,1:NDIM))) then
!                     localParticleCount = localParticleCount + 1
! !                    gr_sbBodyInfo(b) % particles(:,localParticleCount) = SourceBuf(:,k)
!                     !call Grid_updateSolidBodyForces(blkID, b, k, particleposnSB)
!                     exit
!                  endif
!               enddo
!            enddo
!            if (localParticleCount > 0) then
!               send = send + 1
!               sendBufCount = localParticleCount*NPART_PROPS
!               call MPI_ISEND(SourceBuf(1,1), sendBufCount, FLASH_REAL, &
!                    gr_sbBodyInfo(b) % bodyMaster, b, &
!                    gr_meshComm, sendreq(send), ierr)
!            end if
!         endif
!      endif

! !===============================================================
! !================== KPD Change =================================
! !===============================================================
! Call MPI_BARRIER(gr_meshComm, ierr)
! !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! !===============================================================
! !     if (send > 0) then
! !        call MPI_WAITALL(send, sendreq, sstatus, ierr)
! !        if (ierr /= MPI_SUCCESS) then
! !           call Driver_abortFlash("Send MPI_Waitall error")
! !        endif
! !     endif
! !
! !     if (recv > 0) then
! !        call MPI_WAITALL(recv, req, rstatus, ierr)
! !        if (ierr /= MPI_SUCCESS) then
! !           call Driver_abortFlash("Recv MPI_Waitall error")
! !        endif
! !     endif
! !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!   enddo

!   deallocate(sendreq, req, sstatus,rstatus)
  call Timers_stop("body_broadcast_particles")
end subroutine Grid_sbBroadcastParticles
