!!****if* source/Grid/GridStructures/withTriangles/Grid_sbBroadcastParticles
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
!!  * Calculate coordinates of all particles. Calculate the number of particles 
!!  that should be sent from the master to all processors
!!
!!  * Scatter the position of the particles that do not belong to the Master PE
!!    to all processors
!!
!!  * Find which processor the particle belongs to
!!
!! ARGUMENTS 
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine Grid_sbBroadcastParticles()
  use Grid_data, ONLY : gr_meshMe
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, gr_sbPtNumX, &
       gr_sbPtNumY, gr_sbPtNumZ, totalPart
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkBoundBox, &
       Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr
  use Timers_interface, ONLY : Timers_start, Timers_stop
  implicit none
  include "Flash_mpi.h"

  integer :: i, j, k, blkCount, count, numParticles
  integer :: ierr, blkID
  integer :: b, localParticleCount
  real, dimension(2,MDIM) :: boundBox
  real, dimension(MDIM) :: particleposn, particleposnSB
  integer, save :: bufSize
  integer,dimension(MAXBLOCKS) :: listOfBlocks
  real, save, allocatable, dimension(:,:) :: SourceBuf

  call Timers_start("body_broadcast_particles")

!  numParticles = gr_sbPtNumX
!  if (NDIM >= 2) numParticles = numParticles * gr_sbPtNumY
!  if (NDIM == 3) numParticles = numParticles * gr_sbPtNumZ

  allocate(SourceBuf(NPART_PROPS,totalPart))!numParticles)) !SourceBuf contains the particles not in Master

  localParticleCount = 0
  do b = 1, gr_sbNumBodies
     bufSize = 0
     SourceBuf = 0
     if (gr_sbBodyInfo(b) % comm /= MPI_COMM_NULL) then
        ! If master PE, the Master PE gets the positions and the processors
        ! of all particles
        if (gr_sbBodyInfo(b) % myPE == gr_sbBodyInfo(b) % bodyMaster) then
           j = 0
           ! get Blocks in Master PE
           call Grid_getListOfBlocks(LEAF, listOfBlocks, count)
           do i = 1, totalPart!numParticles
              particleposn(IAXIS) = gr_sbBodyInfo(b) % particles(POSX_PART_PROP,i)
              particleposn(JAXIS) = gr_sbBodyInfo(b) % particles(POSY_PART_PROP,i)
              particleposn(KAXIS) = gr_sbBodyInfo(b) % particles(POSZ_PART_PROP,i)
              do blkCount = 1, count
                 blkID = listofBlocks(blkCount)
                 call Grid_getBlkBoundBox(blkID, boundBox)
                 if (all(particleposn(1:NDIM) >= boundBox(LOW,1:NDIM) .and. &
                      particleposn(1:NDIM) < boundBox(HIGH,1:NDIM))) then
                    localParticleCount = localParticleCount + 1
                    call Grid_updateSolidBodyForces(blkID, b, localParticleCount, particleposn)
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
        endif

        ! Broadcast total particle count (totalPart * totalProps)
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

        if (gr_sbBodyInfo(b) % myPE /= gr_sbBodyInfo(b) % bodyMaster) then  !Other processors
           if (bufSize > 0) then
              call Grid_getListOfBlocks(LEAF, listOfBlocks, count) !get blocks                                                                                                                                                               
              do k = 1, j
                 do blkCount = 1, count
                    blkID = listOfBlocks(blkCount)                    
                    call Grid_getBlkBoundBox(blkID, boundBox)
                    particleposnSB(IAXIS) = SourceBuf(POSX_PART_PROP,k)
                    particleposnSB(JAXIS) = SourceBuf(POSY_PART_PROP,k)
                    particleposnSB(KAXIS) = SourceBuf(POSZ_PART_PROP,k)

                    if (all(particleposnSB(1:NDIM) >= boundBox(LOW,1:NDIM) .and. &
                         particleposnSB(1:NDIM) < boundBox(HIGH,1:NDIM))) then
                       localParticleCount = localParticleCount + 1
                       gr_sbBodyInfo(b) % particles(:,localParticleCount) = SourceBuf(:,k)
                       call Grid_updateSolidBodyForces(blkID, b, localParticleCount, particleposnSB)
                       exit
                    endif
                 enddo
              enddo
           endif
        endif
     end if
  enddo
  deallocate(SourceBuf)
  call Timers_stop("body_broadcast_particles")
end subroutine Grid_sbBroadcastParticles
