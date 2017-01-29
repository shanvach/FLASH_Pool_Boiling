!!****if* source/Grid/GridStructures/withTriangles/pointTopoint/bitmap/gr_sbSendForces
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
  use Grid_data, ONLY : gr_meshMe, gr_meshComm, gr_meshNumProcs
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, &
                        gr_sbParticleCount, gr_sbDebug
  use Grid_interface, ONLY : Grid_getBlkBoundBox, Grid_getListOfBlocks
  implicit none
  include "Flash_mpi.h"

  integer :: b, i, recvSize, j, jj, k, p, bufSize, send, recv, &
       particleProc, sendBufCount, recvBufCount, ierr, blkID, gettingFrom
  integer, allocatable, dimension(:) :: sendreq, req
  integer, allocatable, dimension(:,:) :: sstatus, rstatus
  
  integer :: localPart,totalPart

  real, allocatable, dimension(:,:) :: tempparts

  allocate(sstatus(MPI_STATUS_SIZE, gr_meshNumProcs))
  allocate(rstatus(MPI_STATUS_SIZE, gr_meshNumProcs))
  allocate(req(gr_meshNumProcs))
  req = 0
  sstatus = 0
  recvBufCount = NPART_PROPS
  sendBufCount = NPART_PROPS
  send = 0
  recv = 0
  recvSize = 0

  bufSize = count(gr_sbParticleCount(1:gr_sbNumBodies) > 0)
  allocate(sendreq(bufSize))
  sendreq = MPI_REQUEST_NULL

  recv = 0
  do b = 1, gr_sbNumBodies
     if (gr_sbBodyInfo(b) % myPE == gr_sbBodyInfo(b) % bodyMaster) then

        ! Relocate Particles on master in the global numeration GLOB_PART_PROP
        ! Find number of particles in the Master Processor:
        localPart = 0    
        if (associated(gr_sbBodyInfo(b) % particlesPerProc)) then
           do jj=1,size(gr_sbBodyInfo(b)%particlesPerProc,2)
              if(gr_sbBodyInfo(b)%particlesPerProc(1,jj) .eq. gr_sbBodyInfo(b)%bodyMaster) then
                 localPart = gr_sbBodyInfo(b)%particlesPerProc(2,jj)
              endif
           enddo
        endif
        if (gr_sbBodyInfo(b) % sendProcs > 0) then

           if (associated(gr_sbBodyInfo(b) % particlesPerProc)) then
              j = localPart + 1
              recv = 0

              do i = 1, gr_sbBodyInfo(b) % sendProcs !size(gr_sbBodyInfo(b) % particlesPerProc,2)
                 if (gr_sbBodyInfo(b) % particlesPerProc(1,i) /= &
                      gr_sbBodyInfo(b) % bodyMaster) then
                    if (gr_sbBodyInfo(b) % particlesPerProc(2,i) > 0) then
                       recvSize = gr_sbBodyInfo(b) % particlesPerProc(2,i)
                       recvBufCount = recvSize*NPART_PROPS
                       !call MPI_RECV(SourceBuf(1,j), &
                       !     recvBufCount, FLASH_REAL, int(gr_sbBodyInfo(b) % particlesPerProc(1,i)), &
                       !     b, gr_meshComm, rstatus, ierr)
                       recv = recv + 1

                       call MPI_IRECV(gr_sbBodyInfo(b) % particles(1,j), &
                            recvBufCount, FLASH_REAL, int(gr_sbBodyInfo(b) % particlesPerProc(1,i)), &
                            b, gr_meshComm, req(recv), ierr)

!                       print *, "recv", "body", b, "from", gr_sbBodyInfo(b) % particlesPerProc(1,i), &
!                            "count", gr_sbBodyInfo(b) % particlesPerProc(2,i)
!                       do p = 1, recvSize
!                          write(*,'(a,i6,a,3f8.2,a,i6)') "body",b,"Master receiving position", SourceBuf(POSX_PART_PROP:POSZ_PART_PROP,p), &
!                               "from", int(gr_sbBodyInfo(b) % particlesPerProc(1,i))
!                       enddo
                       j = j + recvSize
                    endif
                 endif
              enddo
              !if (gr_sbBodyInfo(b)%sbIsFixed .ne. CONSTANT_ONE) deallocate(gr_sbBodyInfo(b)%particlesPerProc) ! commented by Shizhao to interpolate velocity
           endif
        endif
     endif

     if (gr_sbBodyInfo(b) % myPE /= gr_sbBodyInfo(b) % bodyMaster) then
        gettingFrom = gr_sbParticleCount(b)
        if (gettingFrom > 0) then
           send = 1
           sendBufCount = gettingFrom*NPART_PROPS
!           print *, "send", "body", b, "PE", gr_meshMe, "count", gettingFrom
!           do i = 1, gettingFrom
!              write(*,'(a,i7,a,i7,a,3f8.2)') "body",b, "PE", gr_meshMe,"sending position", &
!                   gr_sbBodyInfo(b) % particles(POSX_PART_PROP:POSZ_PART_PROP,i)
!           enddo
           call MPI_ISEND(gr_sbBodyInfo(b) % particles(1,1), sendBufCount, FLASH_REAL, &
                gr_sbBodyInfo(b) % bodyMaster, b, &
                gr_meshComm, sendreq(send), ierr)
        endif
     endif
     if (send > 0) then
        call MPI_WAITALL(send, sendreq, sstatus, ierr)
        if (ierr /= MPI_SUCCESS) then
           call Driver_abortFlash("Send MPI_Waitall error")
        endif
     endif
     if (recv > 0) then
!        if (recv /= gr_sbBodyInfo(b) % sendProcs) then
!           call Driver_abortFlash("Receive mismatch")
!        end if
        call MPI_WAITALL(recv, req, rstatus, ierr)
        if (ierr /= MPI_SUCCESS) then
           call Driver_abortFlash("Recv MPI_Waitall error")
        endif
     endif
  enddo




  ! Rearrange Particle Positions within Master Processors: -MV 
  do b = 1, gr_sbNumBodies
     if (gr_sbBodyInfo(b) % myPE == gr_sbBodyInfo(b) % bodyMaster) then

        totalPart = gr_sbBodyInfo(b)%totalPart

        if(totalPart .gt. 0) then

           allocate(tempparts(NPART_PROPS,totalPart))

           do jj=1,totalPart

              i = int(gr_sbBodyInfo(b)%particles(GLOB_PART_PROP,jj))

              tempparts(1:NPART_PROPS,i) = gr_sbBodyInfo(b)%particles(1:NPART_PROPS,jj)
          
           enddo

           gr_sbBodyInfo(b)%particles(1:NPART_PROPS,1:totalPart) = tempparts(1:NPART_PROPS,1:totalPart)

           deallocate(tempparts)

        endif
     endif
  enddo

  deallocate(sendreq)
  deallocate(req)
  deallocate(sstatus,rstatus)
!  deallocate(gr_sbParticleCount)
end subroutine gr_sbSendForces
