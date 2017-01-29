!!****if* source/Grid/GridStructures/Grid_updateSolidBodyForces
!!
!! NAME
!!  Grid_updateSolidBodyForces
!!
!! SYNOPSIS
!!
!!  Grid_updateSolidBodyForces()
!!  
!! DESCRIPTION 
!!  
!!  The particle exchange routine
!!
!!  Overview of the algoritm
!!
!!  * Calculate coordinates of all particles. Calculate the number of particles 
!!  that should be sent from the master to all processors
!!
!!  * Scatter the position of the particles that do not belong to the Master PE
!!    to all processors
!!
!!  * Find the block and the processor the particles belongs to
!!
!!  * Force Update: For this unit test we just update the processor ID property.
!!    gr_sbBodyInfo(b) % particles(PROC_PART_PROP,:) = &
!!      real(gr_sbBodyInfo(b) % myPE
!!
!!  * Send the updated particle information to the Master processor. 
!!
!!
!! ARGUMENTS 
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine Grid_updateSolidBodyForces()
  use Grid_data, ONLY : gr_meshMe
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, gr_sbPtNumX, &
       gr_sbPtNumY, gr_sbPtNumZ
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkBoundBox
  implicit none
  include "Flash_mpi.h"

  integer :: b, i, sendSize, recvSize, j, k, blkCount, count, numParticles
  integer :: particleProc, sendBufCount, recvBufCount, ierr, blkID
  real, dimension(2,MDIM) :: boundBox
  integer, save :: bufSize
  real, dimension(MDIM) :: particleposn, particleposnSB
  integer,dimension(MAXBLOCKS) :: listOfBlocks
  integer, allocatable, dimension(:) :: req
  integer, allocatable, dimension(:,:) :: status
  real, allocatable, dimension(:,:) :: SourceBuf

  numParticles = gr_sbPtNumX
  if (NDIM >= 2) numParticles = numParticles * gr_sbPtNumY
  if (NDIM == 3) numParticles = numParticles * gr_sbPtNumZ

  allocate(req(numParticles))
  allocate(SourceBuf(NPART_PROPS,numParticles))
  allocate(status(MPI_STATUS_SIZE, numParticles))

  status = 0

  do b = 1, gr_sbNumBodies
     bufSize = 0
     req(:)=0
     SourceBuf = 0
     if (gr_sbBodyInfo(b) % comm /= MPI_COMM_NULL) then
        ! If master PE, the Master PE gets the positions and the processors
        ! of all particles
        if (gr_sbBodyInfo(b) % myPE == gr_sbBodyInfo(b) % bodyMaster) then
           j = 0
           ! get Blocks in Master PE
           call Grid_getListOfBlocks(LEAF, listOfBlocks, count)
           do i = 1, numParticles
              particleposn(IAXIS) = gr_sbBodyInfo(b) % particles(POSX_PART_PROP,i)
              particleposn(JAXIS) = gr_sbBodyInfo(b) % particles(POSY_PART_PROP,i)
              particleposn(KAXIS) = gr_sbBodyInfo(b) % particles(POSZ_PART_PROP,i)
              do blkCount = 1, count
                 blkID = listofBlocks(blkCount)
                 call Grid_getBlkBoundBox(blkID, boundBox)
                 if (all(particleposn(1:NDIM) >= boundBox(LOW,1:NDIM) .and. &
                      particleposn(1:NDIM) < boundBox(HIGH,1:NDIM))) then
                    !particle in Master PE
                    exit
                 endif
              enddo
              ! If not in master PE, no. of particles to be sent by Master PE
              if (blkCount > count) then
                 bufSize = bufSize + NPART_PROPS 
                 j = j + 1
                 SourceBuf(1:NPART_PROPS,j) = &
                      gr_sbBodyInfo(b) % particles(1:NPART_PROPS,i)
              endif
           enddo

           recvSize = 0
           !if particles sent by master to other processors, then master
           !needs to receive it
           if (j > 0) then
              do k = 1, numParticles
                 recvBufCount = NPART_PROPS
                 particleProc = gr_sbBodyInfo(b) % particles(PROC_PART_PROP,k)
                 ! particle does not belong to master processor
                 if (particleProc /= gr_sbBodyInfo(b) % bodyMaster) then
                    recvSize = recvSize + 1
                    call MPI_IRECV(gr_sbBodyInfo(b) % particles(1,k), &
                         recvBufCount, FLASH_REAL, MPI_ANY_SOURCE, 1, &
                         gr_sbBodyInfo(b) % comm, req(recvSize), ierr)
                 endif
              enddo
           endif
        endif

        ! Broadcast total particle count (numParticles * totalProps)
        call MPI_BCAST(bufSize, 1, FLASH_INTEGER, &
             gr_sbBodyInfo(b) % bodyMaster, gr_sbBodyInfo(b) % comm, ierr)
        if (bufSize > 0) then
           ! Broadcast particle positions to all processors
           call MPI_BCAST(SourceBuf, bufSize, FLASH_REAL, &
                gr_sbBodyInfo(b) % bodyMaster, gr_sbBodyInfo(b) % comm, ierr)
           ! Broadcast particle count (no of particles)
           call MPI_BCAST(j, 1, FLASH_INTEGER, gr_sbBodyInfo(b) % bodyMaster, &
                gr_sbBodyInfo(b) % comm, ierr)
        endif

        ! If not Master PE, find the block and the PE the particle belongs to
        if (gr_sbBodyInfo(b) % myPE /= gr_sbBodyInfo(b) % bodyMaster) then
           if (bufSize > 0) then
              call Grid_getListOfBlocks(LEAF, listOfBlocks, count) !get blocks
              do k = 1, j
                 do blkCount = 1, count
                    blkID = listOfBlocks(blkCount)
                    call Grid_getBlkBoundBox(blkID, boundBox)
                    sendSize = 0
                    particleposnSB(IAXIS) = SourceBuf(POSX_PART_PROP,k)
                    particleposnSB(JAXIS) = SourceBuf(POSY_PART_PROP,k)
                    particleposnSB(KAXIS) = SourceBuf(POSZ_PART_PROP,k)

                    if (all(particleposnSB(1:NDIM) >= boundBox(LOW,1:NDIM) .and. &
                         particleposnSB(1:NDIM) < boundBox(HIGH,1:NDIM))) then
                       !update processor ID.
                       SourceBuf(PROC_PART_PROP,k) = real(gr_sbBodyInfo(b) % myPE)
                       sendBufCount = NPART_PROPS
                       sendSize = sendSize + 1
                       call MPI_ISEND(SourceBuf(1,k), sendBufCount, FLASH_REAL, &
                            gr_sbBodyInfo(b) % bodyMaster, 1, &
                            gr_sbBodyInfo(b) % comm, req(sendSize), ierr)
                       exit
                    endif
                 enddo
                 ! If particles have been sent to Master PE, block until all
                 ! sends are complete
                 if (sendSize > 0) then
                    call MPI_WAITALL(sendSize, req, status, ierr)
                 end if
              end do
           endif
        else  ! If Master PE
           ! If Master PE needs to receive messages, block until all
           ! receives complete
           if (recvSize > 0) then
              call MPI_WAITALL(recvSize, req, status, ierr) 
           end if
        endif
     end if
  enddo
end subroutine Grid_updateSolidBodyForces
