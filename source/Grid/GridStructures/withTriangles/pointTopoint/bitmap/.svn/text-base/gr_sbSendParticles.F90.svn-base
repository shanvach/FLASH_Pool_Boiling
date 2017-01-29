!!****if* source/Grid/GridStructures/withTriangles/pointTopoint/bitmap/gr_sbSendParticles
!!
!! NAME
!!  gr_sbSendParticles
!!
!! SYNOPSIS
!!  gr_sbSendParticles()
!!
!! DESCRIPTION
!!
!!  * Actually send the particles.
!!  * The subroutine gr_sbSendParticleCount has already informed
!!    the slave processors how many particles they will receive.
!!  * Sorts the particles according to the destination processor.
!!  * Sends particles to appropriate slave processor
!!  * Slave Processors check if needs to receive particles. If yes, then receives 
!!    particles from master
!!
!! ARGUMENTS
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine gr_sbSendParticles
  use Grid_data, ONLY : gr_meshComm, gr_meshNumProcs, gr_meshMe
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, totalPart, sumPart, gr_sbDebug, &
                        gr_sbParticleCount, solid_body, gr_sbFirstCall, gr_sbStencil
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use ut_sortInterface, ONLY : ut_sortOnProcs
  use Logfile_interface, ONLY : Logfile_open, Logfile_close
  implicit none
  include "Flash_mpi.h"
  type(solid_body), pointer :: bodyInfo
  real, dimension(NPART_PROPS,sumPart) :: TempBuf
  integer :: i, j, l, b, countProc, gettingFrom, p, &
       ierr, sendCount, recvCount, recvBufCount, sendBufCount, &
       iSendReq, iRecvReq, logUnit, totalSendMsg, totalRecvMsg
  integer, dimension(0:gr_meshNumProcs-1) :: perProc, logProc
  logical, parameter :: localLogFile = .true.
  integer, allocatable, dimension(:) :: sreq, rreq
  integer, allocatable, dimension(:,:) :: sstatus, rstatus

  call Timers_start("particle_transfer")
  iSendReq = 0
  iRecvReq = 0

  !Count the total number of messages that we will send so we can 
  !allocate the send status and request objects appropriately.
  totalSendMsg = 0
  do b = 1, gr_sbNumBodies
     bodyInfo => gr_sbBodyInfo(b)
     if (bodyInfo % myPE == bodyInfo % bodyMaster) then
        if (associated(bodyInfo % particlesPerProc)) then
           do j = 1, size(bodyInfo % particlesPerProc,2)
              if (bodyInfo % particlesPerProc(1,j) /= bodyInfo % bodyMaster) then
                 totalSendMsg = totalSendMsg + 1
              end if
           end do
!           deallocate(bodyInfo % particlesPerProc) !No longer needed.
        end if
     end if
  end do
  allocate(sstatus(MPI_STATUS_SIZE,totalSendMsg))
  allocate(sreq(totalSendMsg))
  sreq = MPI_REQUEST_NULL

  !Count the total number of messages that we will receive so we can
  !allocate the receive status and request objects appropriately.
  totalRecvMsg = count(gr_sbParticleCount(1:gr_sbNumBodies) > 0)
  allocate(rstatus(MPI_STATUS_SIZE,totalRecvMsg))
  allocate(rreq(totalRecvMsg))
  rreq = MPI_REQUEST_NULL

  !The master process sorts the particles array and then sends
  !all off-processor particles
  do b = 1, gr_sbNumBodies

     ! IF fixed body and not the first call:
     if ((gr_sbBodyInfo(b)%sbIsFixed .eq. CONSTANT_ONE) .and. (gr_sbFirstCall .eq. CONSTANT_ZERO)) cycle

     bodyInfo => gr_sbBodyInfo(b)
     if (bodyInfo % myPE == bodyInfo % bodyMaster) then
        totalPart = bodyInfo%totalPart
        !perProc: Array of size = no of processors. Each entry is the count of particles
        !         destined for the corresp processor
        !logProc: Array of size = no of processors. Value = 1 if perProc > 0
        !countProc: count of processors that need to recieve particles
        call ut_sortOnProcs(totalPart, NPART_PROPS, PROC_PART_PROP, &
             gr_meshNumProcs, bodyInfo % particles, TempBuf, &
             perProc, logProc, countProc)       
        bodyInfo % sendProcs = countProc
        if(countProc>0) then
           j=0
           l=1
           do i=1,countProc
              do while(perProc(j)==0)
                 j=j+1
              end do
              sendCount = perProc(j)
              sendBufCount = sendCount*NPART_PROPS
              
              if (j /= bodyInfo % bodyMaster) then
                 if (gr_sbDebug) then
                    write(logUnit,*) "Body", b, "master sending", &
                         sendBufCount, " to ", j
                 end if

                 iSendReq = iSendReq + 1
                 if (iSendReq > totalSendMsg) then
                    call Driver_abortFlash("Send counter overflow")
                 end if
!                 do k = l, sendCount
!                    write(*,'(a,i7,a,i7,a,3f8.2,a,i6)') "body",b, "master", gr_sbBodyInfo(b) % bodyMaster,"sending position", &
!                          bodyInfo % particles(POSX_PART_PROP:POSZ_PART_PROP,k), "to", j
!                enddo
                 call MPI_Isend(bodyInfo % particles(1,l),sendBufCount,&
                      FLASH_REAL,j,b,gr_meshComm,sreq(iSendReq),ierr)
                 if (ierr /= MPI_SUCCESS) then
                    call Driver_abortFlash("MPI_Isend error")
                 endif
              end if
              j = j+1
              l = l+sendCount
           end do
        end if
     else
        gettingFrom = gr_sbParticleCount(b)
        bodyInfo % totalPart = 0
        if (gettingFrom > 0) then
           bodyInfo % totalPart = gettingFrom
           allocate(bodyInfo % particles(NPART_PROPS,gettingFrom))
           if (bodyInfo%sbIsFixed .eq. CONSTANT_ONE) then
             allocate(bodyInfo%ielem(gr_sbStencil,MDIM,MDIM,gettingFrom))
             allocate(bodyInfo%phile(gr_sbStencil,NDIM+1,MDIM,gettingFrom))
             bodyInfo%ielem = NONEXISTENT
             bodyInfo%phile = 0.
           endif
           recvCount = gettingFrom
           recvBufCount = recvCount * NPART_PROPS
           if (gr_sbDebug) then
              write(logUnit,*) "Body", b, "slave receiving", recvBufCount, &
                   " from ", gr_sbBodyInfo(b) % bodyMaster
           end if

           iRecvReq = iRecvReq + 1
           if (iRecvReq > totalRecvMsg) then
              call Driver_abortFlash("Receive counter overflow")
           end if
           call MPI_Irecv(bodyInfo % particles(1,1),recvBufCount,&
                FLASH_REAL,bodyInfo % bodyMaster,b,gr_meshComm,&
                rreq(iRecvReq),ierr)
            if (ierr /= MPI_SUCCESS) then
              call Driver_abortFlash("MPI_Irecv error")
           endif
        end if
     end if

  end do

  if (iSendReq > 0) then
     call MPI_Waitall(iSendReq, sreq, sstatus, ierr)
     if (ierr /= MPI_SUCCESS) then
        call Driver_abortFlash("Send MPI_Waitall error")
     endif
  end if
  if (iRecvReq > 0) then
     call MPI_Waitall(iRecvReq, rreq, rstatus, ierr)
     do b = 1, gr_sbNumBodies
        gettingFrom = gr_sbParticleCount(b)
        do p = 1, gettingFrom
           if (gr_sbBodyInfo(b) % particles(PROC_PART_PROP,p) /= gr_meshMe) then
              print *, "ERROR"
           endif
        end do
     end do
!     do b = 1, gr_sbNumBodies
!        gettingFrom = gr_sbParticleCount(b)
!        do p = 1, gettingFrom
!           write(*,'(a, i6, a,3f8.2,a,i6)') "myID", gr_meshMe, "receiving position", gr_sbBodyInfo(b) % particles(POSX_PART_PROP:POSZ_PART_PROP,p), &
!                "from", gr_sbBodyInfo(b) % bodyMaster
!        enddo        
!     enddo
     if (ierr /= MPI_SUCCESS) then
        call Driver_abortFlash("Recv MPI_Waitall error")
     endif
  end if

  deallocate(sstatus, rstatus, sreq, rreq)

  call Timers_stop("particle_transfer")
end subroutine gr_sbSendParticles
