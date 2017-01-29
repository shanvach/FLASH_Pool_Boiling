!!****if* source/Simulation/SimulationMain/INavierStokes/2D/IB_Cyl_parallelIBVP/old/gr_sbSendForces
!!
!! NAME
!!  gr_sbSendForces
!!
!! SYNOPSIS
!!
!!  gr_sbSendForces()
!!  
!! DESCRIPTION 
!!  
!!  The particle exchange routine
!!
!!  Overview of the algoritm
!!
!!  The slave processors send the updated particle information 
!!  to the Master processor. 
!!
!!
!! ARGUMENTS 
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine gr_sbSendForces()
  use Grid_data, ONLY : gr_meshMe, gr_meshComm
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, &
       totalPart, gr_sbParticleCount, gr_sbDebug
  use Grid_interface, ONLY : Grid_getBlkBoundBox, Grid_getListOfBlocks
  implicit none
  include "Flash_mpi.h"

  integer :: b, i, sendSize, recvSize, j, k, count, &
       particleProc, sendBufCount, recvBufCount, ierr, blkID, gettingFrom
  integer, allocatable, dimension(:) :: req
  integer, allocatable, dimension(:,:) :: status
  real, dimension(2,MDIM) :: boundBox
  integer, dimension(MAXBLOCKS) :: listOfBlocks
  real, dimension(MDIM) :: pos
  

  recvBufCount = NPART_PROPS
  sendBufCount = NPART_PROPS

  do b = 1, gr_sbNumBodies

     totalPart = gr_sbBodyInfo(b) % totalPart
 
     allocate(req(totalPart))
     allocate(status(MPI_STATUS_SIZE, totalPart))
     status = 0
     req(:)=0

     if (gr_sbBodyInfo(b) % myPE == gr_sbBodyInfo(b) % bodyMaster) then
        recvSize = 0
        !if particles sent by master to other processors, then master
        !needs to receive it
!        print *, "body", b, "procs", gr_sbBodyInfo(b) % sendProcs
        if (gr_sbBodyInfo(b) % sendProcs > 0) then
           do k = 1, totalPart
              particleProc = gr_sbBodyInfo(b) % particles(PROC_PART_PROP,k)
              ! particle does not belong to master processor
              if (particleProc /= gr_sbBodyInfo(b) % bodyMaster) then
                 recvSize = recvSize + 1
                 call MPI_IRECV(gr_sbBodyInfo(b) % particles(1,k), &
                      recvBufCount, FLASH_REAL, MPI_ANY_SOURCE, b, &
                      gr_meshComm, req(recvSize), ierr)
              endif
           enddo
        endif
     endif

     if (gr_sbBodyInfo(b) % myPE /= gr_sbBodyInfo(b) % bodyMaster) then
        gettingFrom = gr_sbParticleCount(b)
        if (gettingFrom > 0) then
           sendSize = 0
           do k = 1, gettingFrom
              sendSize = sendSize + 1
              call MPI_ISEND(gr_sbBodyInfo(b) % particles(1,k), sendBufCount, FLASH_REAL, &
                   gr_sbBodyInfo(b) % bodyMaster, b, &
                   gr_meshComm, req(sendSize), ierr)
           enddo
           ! If particles have been sent to Master PE, block until all
           ! sends are complete
           if (sendSize > 0) then
              call MPI_WAITALL(sendSize, req, status, ierr)
           end if
           deallocate(gr_sbBodyInfo(b) % particles)
        endif
     else  ! If Master PE
        ! If Master PE needs to receive messages, block until all
        ! receives complete
        if (recvSize > 0) then
           call MPI_WAITALL(recvSize, req, status, ierr)            
        end if
        call Grid_getListOfBlocks(LEAF, listOfBlocks, count)

        !Check sends and receives
        if (gr_sbDebug) then
           do k = 1, totalPart
              particleProc = gr_sbBodyInfo(b) % particles(PROC_PART_PROP,k)
              if (particleProc /= gr_sbBodyInfo(b) % bodyMaster) then
                 pos(IAXIS) = gr_sbBodyInfo(b) % particles(POSX_PART_PROP,k)
                 pos(JAXIS) = gr_sbBodyInfo(b) % particles(POSY_PART_PROP,k)
                 pos(KAXIS) = gr_sbBodyInfo(b) % particles(POSZ_PART_PROP,k)
                 do i = 1, count
                    blkID = listOfBlocks(i)
                    call Grid_getBlkBoundBox(blkID, boundBox)
                    if (all(&
                         pos(1:NDIM) >= boundBox(LOW,1:NDIM) .and. &
                         pos(1:NDIM) < boundBox(HIGH,1:NDIM))) then
                       print *, "error in sends and receives"
                    endif
                 enddo
              endif
           enddo
        endif
     endif

     deallocate(req)
     deallocate(status)


  enddo
  deallocate(gr_sbParticleCount)
end subroutine gr_sbSendForces
