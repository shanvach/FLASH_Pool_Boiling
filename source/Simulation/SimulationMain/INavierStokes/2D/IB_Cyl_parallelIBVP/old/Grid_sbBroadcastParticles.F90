!!****if* source/Simulation/SimulationMain/INavierStokes/2D/IB_Cyl_parallelIBVP/old/Grid_sbBroadcastParticles
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
!! ARGUMENTS 
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine Grid_sbBroadcastParticles()
  use Grid_data, ONLY : gr_meshMe, gr_meshComm
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, &
      totalPart, sumPart
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkBoundBox, &
       Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr
  use Timers_interface, ONLY : Timers_start, Timers_stop
  implicit none
  include "Flash_mpi.h"

  integer :: i, j, k, blkCount, count
  integer :: ierr, blkID
  integer :: b, localParticleCount
  integer,save :: particleProc
  real, dimension(2,MDIM) :: boundBox
  real, dimension(MDIM) :: particleposn, particleposnSB
  integer:: bufSize
  integer,dimension(MAXBLOCKS) :: listOfBlocks
  real, dimension(NPART_PROPS,sumPart) :: SourceBuf
  real :: particleData(NPART_PROPS)

  call Timers_start("body_broadcast_particles")

  particleposn = 0.

  do b = 1, gr_sbNumBodies
     bufSize = 0
     SourceBuf = 0
        ! If master PE, the Master PE gets the positions and the processors
        ! of all particles
     if (gr_sbBodyInfo(b) % myPE == gr_sbBodyInfo(b) % bodyMaster) then
        localParticleCount = 0
        j = 0
           ! get Blocks in Master PE
        call Grid_getListOfBlocks(LEAF, listOfBlocks, count)

        totalPart = gr_sbBodyInfo(b) % totalPart

        do i = 1, totalPart
           particleposn(IAXIS) = gr_sbBodyInfo(b) % particles(POSX_PART_PROP,i)
           if(NDIM > 1) then
              particleposn(JAXIS) = gr_sbBodyInfo(b) % particles(POSY_PART_PROP,i)
           endif
           if(NDIM > 2) then
              particleposn(KAXIS) = gr_sbBodyInfo(b) % particles(POSZ_PART_PROP,i)
           endif
           do blkCount = 1, count
              blkID = listofBlocks(blkCount)
              call Grid_getBlkBoundBox(blkID, boundBox)
              if (all(particleposn(1:NDIM) >= boundBox(LOW,1:NDIM) .and. &
                   particleposn(1:NDIM) < boundBox(HIGH,1:NDIM))) then
!                 localParticleCount = localParticleCount + 1
                 gr_sbBodyInfo(b) % particles(BLK_PART_PROP,i) = blkID 
                 gr_sbBodyInfo(b) % particles(PROC_PART_PROP,i) = gr_sbBodyInfo(b) % myPE
                 particleData = gr_sbBodyInfo(b) % particles(1:NPART_PROPS,i)



!                 call Grid_updateSolidBodyForces(blkID, particleData)
                 !call Grid_updateSolidBodyForces(blkID, b, i, particleposn)
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
     
        ! Broadcast total particle count (numParticles * totalProps)
     call MPI_BCAST(bufSize, 1, FLASH_INTEGER, &
          gr_sbBodyInfo(b) % bodyMaster, MPI_COMM_WORLD, ierr)
     if (bufSize > 0) then
        ! Broadcast particle positions to all processors
        call MPI_BCAST(SourceBuf, bufSize, FLASH_REAL, &
             gr_sbBodyInfo(b) % bodyMaster, MPI_COMM_WORLD, ierr)
           ! Broadcast particle count (no of particles)
        call MPI_BCAST(j, 1, FLASH_INTEGER, gr_sbBodyInfo(b) % bodyMaster, &
             MPI_COMM_WORLD, ierr)
     endif

     if (gr_sbBodyInfo(b) % myPE /= gr_sbBodyInfo(b) % bodyMaster) then  !Other processors
        localParticleCount = 0
        if (bufSize > 0) then
           call Grid_getListOfBlocks(LEAF, listOfBlocks, count) !get blocks  
           do k = 1, j
              do blkCount = 1, count
                 blkID = listOfBlocks(blkCount)                    
                 call Grid_getBlkBoundBox(blkID, boundBox)
                 particleposnSB(IAXIS) = SourceBuf(POSX_PART_PROP,k)
                 if(NDIM > 1) then
                    particleposnSB(JAXIS) = SourceBuf(POSY_PART_PROP,k)
                 endif
                 if(NDIM > 2) then
                    particleposnSB(KAXIS) = SourceBuf(POSZ_PART_PROP,k)
                 endif
                 
                 if (all(particleposnSB(1:NDIM) >= boundBox(LOW,1:NDIM) .and. &
                      particleposnSB(1:NDIM) < boundBox(HIGH,1:NDIM))) then
!                    localParticleCount = localParticleCount + 1
!                    gr_sbBodyInfo(b) % particles(:,localParticleCount) = SourceBuf(:,k)
                    particleData = Sourcebuf(1:NPART_PROPS,k)
!                    call Grid_updateSolidBodyForces(blkID, particleData)
                    !call Grid_updateSolidBodyForces(blkID, b, k, particleposnSB)
                    exit
                 endif
              enddo
           enddo
        endif
     endif
  enddo
  call Timers_stop("body_broadcast_particles")
end subroutine Grid_sbBroadcastParticles
