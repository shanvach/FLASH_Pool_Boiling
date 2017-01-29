!!****if* source/Grid/GridStructures/OneSideComm/Grid_updateSolidBodyForces
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
!!  The particle exchange routine (MPI-2 comm)
!!
!!  Overview of the algoritm
!!
!!  * Master processor creates a local window and buffer that will have the
!!  update processor ID of the particle
!!
!!  * Find the block and the processor the particles belongs to
!!
!!  * Force Update: For this unit test we just update the processor ID property. 
!! gr_sbBodyInfo(b) % particles(PROC_PART_PROP,:) = real(gr_sbBodyInfo(b) % myPE
!!
!!  * Send the updated particle information to the Master processor using
!!    one-sided communication
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
  use Timers_interface, ONLY : Timers_start, Timers_stop
  implicit none
  include "Flash_mpi.h"

  integer :: b, k, blkCount, count, numParticles, particleProc, ierr, blkID, win
  real, dimension(2,MDIM) :: boundBox
  real, dimension(MDIM) :: particleposn
  integer,dimension(MAXBLOCKS) :: listOfBlocks
  real, allocatable, dimension(:,:) :: SourceBuf
  real, allocatable, dimension(:) :: ProcID
  integer (kind=MPI_ADDRESS_KIND) :: sizeofreal, lowerbound, targetDisp

  call Timers_start("body_forces")
  call MPI_Type_get_extent(FLASH_REAL, lowerbound, sizeofreal, ierr)

  numParticles = gr_sbPtNumX
  if (NDIM >= 2) numParticles = numParticles * gr_sbPtNumY
  if (NDIM == 3) numParticles = numParticles * gr_sbPtNumZ

  allocate(SourceBuf(NPART_PROPS,numParticles+1))
  !local buffer of Master PE that gets the updated processor ID of particles
  allocate(ProcID(numParticles))

  ProcID = NONEXISTENT

  do b = 1, gr_sbNumBodies  
!     print *, "body no", b
     if (gr_sbBodyInfo(b) % comm /= MPI_COMM_NULL) then

        call MPI_Bcast(gr_sbBodyInfo(b) % particles, numParticles*NPART_PROPS, &
             FLASH_REAL, gr_sbBodyInfo(b) % bodyMaster, gr_sbBodyInfo(b) % comm, ierr)

       if (gr_sbBodyInfo(b) % myPE == gr_sbBodyInfo(b) % bodyMaster) then
          do k = 1, numParticles
             particleProc = gr_sbBodyInfo(b) % particles(PROC_PART_PROP,k)
             if (particleProc == gr_sbBodyInfo(b) % bodyMaster) then
                !If particle is in the Master PE, no change in Processor ID
                ProcID(k) = gr_sbBodyInfo(b) % bodyMaster
                !print *, "Master receiveing particle", k, "position", &
                !gr_sbBodyInfo(b) % particles(POSX_PART_PROP:POSY_PART_PROP, k), &
                !"PE", gr_sbBodyInfo(b) % particles(PROC_PART_PROP, k)
             endif
          enddo
          !Master PE creates a local window (memory) that will get the
          !updated Proc ID
          call MPI_Win_create(ProcID, numParticles*sizeofreal, sizeofreal, &
               MPI_INFO_NULL, gr_sbBodyInfo(b) % comm, win, ierr) 
       else
          !Other processors create Null memory. Since put operation in use,
          !data sent need not be in the local window 
          call MPI_Win_create(MPI_BOTTOM, 0, sizeofreal, MPI_INFO_NULL, &
               gr_sbBodyInfo(b) % comm, win, ierr) 
       endif
       
       !Prepare for 1st put
       call MPI_Win_fence(0, win, ierr)

       ! If not Master PE, find the block and the PE the particle belongs to
       if (gr_sbBodyInfo(b) % myPE /= gr_sbBodyInfo(b) % bodyMaster) then
          call Grid_getListOfBlocks(LEAF, listOfBlocks, count) !get blocks
          do k = 1, numParticles
             SourceBuf(1:NPART_PROPS,k) = &
                  gr_sbBodyInfo(b) % particles(1:NPART_PROPS,k)
             particleposn(IAXIS) = gr_sbBodyInfo(b) % particles(POSX_PART_PROP,k)
             particleposn(JAXIS) = gr_sbBodyInfo(b) % particles(POSY_PART_PROP,k)
             particleposn(KAXIS) = gr_sbBodyInfo(b) % particles(POSZ_PART_PROP,k)
             do blkCount = 1, count
                blkID = listOfBlocks(blkCount)
                call Grid_getBlkBoundBox(blkID, boundBox)

                if (all(particleposn(1:NDIM) >= boundBox(LOW,1:NDIM) .and. &
                     particleposn(1:NDIM) < boundBox(HIGH,1:NDIM))) then
                   !update processor ID.
                   SourceBuf(PROC_PART_PROP,k) = real(gr_sbBodyInfo(b) % myPE)
                   targetDisp = k-1
                   call MPI_Put(SourceBuf(PROC_PART_PROP, k), 1, FLASH_REAL, &
                        gr_sbBodyInfo(b) % bodyMaster, targetDisp, 1, &
                        FLASH_REAL, win, ierr) ! write to Master PE memory
                   !print *, "Slave processor", gr_sbBodyInfo(b) % myPE, &
                   !     " sending particle processor", SourceBuf(PROC_PART_PROP,k), &
                   !     "to position", int(k-1)
                   exit ! Exit when particle found in a block. Go to next particle
                endif
             enddo
          end do
       endif

       !All puts now done
       call MPI_Win_fence(0, win, ierr)

       ! If master PE, the Master PE gets the positions and the
       ! processors of all particles
       if (gr_sbBodyInfo(b) % myPE == gr_sbBodyInfo(b) % bodyMaster) then
          do k = 1, numParticles
             !print *, int(ProcID(k))
             particleProc = gr_sbBodyInfo(b) % particles(PROC_PART_PROP,k)
             if (particleProc /= gr_sbBodyInfo(b) % bodyMaster) then
                gr_sbBodyInfo(b) % particles(PROC_PART_PROP,k) = real(ProcID(k))
             endif
          end do
       end if
       call MPI_Win_free(win, ierr)
    end if
 enddo
 deallocate(SourceBuf)
 deallocate(ProcID)

 call Timers_stop("body_forces")
end subroutine Grid_updateSolidBodyForces
