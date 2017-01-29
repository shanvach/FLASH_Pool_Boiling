!!****if* source/Grid/GridStructures/withTriangles/pointTopoint/bitmap/gr_sbSendPosn
!!
!! NAME
!!  gr_sbSendPosn
!!
!! SYNOPSIS
!!
!!  gr_sbSendPosn()
!!  
!! DESCRIPTION 
!!  
!!  The particle broadcast routine
!!
!!  Overview of the algoritm
!!
!!  * Find if partice belongs to Master PE. If not, store it in a buffer.
!!  * Sorts the particles according to the destination processor.
!!  * Send the no. of particles each slave processor needs to receive from 
!!    master to each appropriate slave processor.   
!!  * Sends particles to appropriate slave processor
!!  * Slave Processors check if needs to receive particles. If yes, then receives 
!!    particles from master.
!!
!! ARGUMENTS 
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine gr_sbSendPosn()
  use Timers_interface, ONLY : Timers_start, Timers_stop
  implicit none
!  call Timers_start("send_particles_posn")
  call gr_sbStoreParticlesPerProc()
  call gr_sbSendParticleCount()
  call gr_sbSendParticles()
!  call Timers_start("gr_sbUpDateForces")
!  call gr_sbUpdateForces()
!  call gr_sbUpdateForces_Lag()
!  call gr_sbUpdateForces_loc()
!  call Timers_stop("gr_sbUpDateForces")
!  call Timers_stop("send_particles_posn")
end subroutine gr_sbSendPosn


!This subroutine steps through the particles array and
!counts how many particles are destined for each processor.
!The processor ID and particle count is stored in
!an array in the gr_sbBodyInfo data structure.
!
!Side effects:
!Modifies the gr_sbBodyInfo data structure.
!Allocates and initializes the gr_sbParticleCount array.
!subroutine gr_sbStoreParticlesPerProc
!  use Timers_interface, ONLY : Timers_start, Timers_stop
!  use Grid_data, ONLY : gr_meshMe, gr_meshComm, gr_meshNumProcs
!  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, gr_sbDebug, &
!       totalPart, solid_body, gr_sbParticleCount
!  implicit none
!  integer, dimension(0:gr_meshNumProcs-1) :: perProc
!  type(solid_body), pointer :: bodyInfo
!  integer :: i, j, b, procID, proc, countProcs

!  call Timers_start("store_particles_per_proc")
!  do b = 1, gr_sbNumBodies
!     bodyInfo => gr_sbBodyInfo(b)
!     if (bodyInfo % myPE == bodyInfo % bodyMaster) then
!        perProc = 0
!        countProcs = 0
!        do i = 1, totalPart
!           procID = int(bodyInfo % particles(PROC_PART_PROP,i))
!           if (procID >= 0) then
!              if (perProc(procID) == 0) then
!                 countProcs = countProcs + 1 !counts the no of processors that need to receive particles
!              end if
!              perProc(procID) = perProc(procID) + 1
!           end if
!        end do

!        if (countProcs > 0) then
!           allocate(bodyInfo % particlesPerProc(1:2,1:countProcs))
!           proc = 0
!           do j = 0, gr_meshNumProcs-1
!              if (perProc(j) > 0) then
!                 proc = proc + 1
!                 bodyInfo % particlesPerProc(1,proc) = j
!                 bodyInfo % particlesPerProc(2,proc) = perProc(j)
!              end if
!          end do
!        else
!           nullify(bodyInfo % particlesPerProc)
!        end if
!     end if
! end do

  !We allocate and initialize in this subroutine and not
  !gr_sbSendParticleCount to ensure the compiler flushes
  !local initial values to memory before we start remote writes.
!  allocate(gr_sbParticleCount(1:gr_sbNumBodies))
!  gr_sbParticleCount = 0
!  call Timers_stop("store_particles_per_proc")
!end subroutine gr_sbStoreParticlesPerProc


!We use the processor ID and particle count information
!collected in gr_sbStoreParticlesPerProc to inform
!the required slave processors how many particles they will
!receive.  The communication depends upon MPI-2 and is
!completely contained between two synchronization fences.
!
!Side effects:
!Updates an array named gr_sbParticleCount
!subroutine gr_sbSendParticleCount
!  use Timers_interface, ONLY : Timers_start, Timers_stop
!  use Grid_data, ONLY : gr_meshMe, gr_meshComm, gr_meshNumProcs
!  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, gr_sbDebug, &
!       totalPart, solid_body, gr_sbParticleCount
!  implicit none
!  include "Flash_mpi.h"
!  type(solid_body), pointer :: bodyInfo
!  integer :: win, ierr, b, j, targetRank, dispUnit
!  integer (kind=MPI_ADDRESS_KIND) :: lowerbound, intSize, winSize, targetDisp

!  call Timers_start("send_particle_count")
!  call MPI_Type_get_extent(MPI_INTEGER, lowerbound, intSize, ierr)
!  winSize = intSize * gr_sbNumBodies
!  dispUnit = intSize
!  call MPI_Win_create(gr_sbParticleCount, winSize, dispUnit, &
!       MPI_INFO_NULL, gr_meshComm, win, ierr)

!  call MPI_Win_fence(0, win, ierr)
!  do b = 1, gr_sbNumBodies
!     bodyInfo => gr_sbBodyInfo(b)
!     if (bodyInfo % myPE == bodyInfo % bodyMaster) then
!        if (associated(bodyInfo % particlesPerProc)) then
!           do j = 1, size(bodyInfo % particlesPerProc,2)
!              targetRank = bodyInfo % particlesPerProc(1,j)
!              if (targetRank /= bodyInfo % bodyMaster) then
!                 !print *, "PE", bodyInfo % bodyMaster, "sending", &
!                 !     bodyInfo % particlesPerProc(2,j), "to", &
!                 !     targetRank, "for body:", b
!                 targetDisp = b-1
                 !Note that MPI_Put behaves like MPI_Isend and so
                 !communication may be deferred until the MPI_Win_fence.
!                 call MPI_Put(bodyInfo % particlesPerProc(2,j), 1, MPI_INTEGER, &
!                      targetRank, targetDisp, 1, MPI_INTEGER, win, ierr)
!              end if
!           end do
!        end if
!     end if
!  end do
!  call MPI_Win_fence(0, win, ierr)

!  call MPI_Win_free(win, ierr)
!  call Timers_stop("send_particle_count")
!end subroutine gr_sbSendParticleCount


!Actually send the particles.
!The subroutine gr_sbSendParticleCount has already informed
!the slave processors how many particles they will receive.
!subroutine gr_sbSendParticles
!  use Grid_data, ONLY : gr_meshComm, gr_meshNumProcs
!  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, totalPart, gr_sbDebug, &
!       gr_sbParticleCount, solid_body
!  use Timers_interface, ONLY : Timers_start, Timers_stop
!  use ut_sortInterface, ONLY : ut_sortOnProcs
!  use Logfile_interface, ONLY : Logfile_open, Logfile_close
!  implicit none
!  include "Flash_mpi.h"
!  type(solid_body), pointer :: bodyInfo
!  real, dimension(NPART_PROPS,totalPart) :: TempBuf
!  integer :: i, j, l, b, countProc, totalProcPart, gettingFrom, countAndProc, &
!       ierr, sendCount, recvCount, recvBufCount, sendBufCount, &
!       iSendReq, iRecvReq, logUnit, totalSendMsg, totalRecvMsg
!  integer, dimension(0:gr_meshNumProcs-1) :: perProc, logProc, fromProcs
!  logical, parameter :: localLogFile = .true.
!  integer, allocatable, dimension(:) :: sreq, rreq
!  integer, allocatable, dimension(:,:) :: sstatus, rstatus

!  call Timers_start("particle_transfer")
!  iSendReq = 0
!  iRecvReq = 0

  !Count the total number of messages that we will send so we can
  !allocate the send status and request objects appropriately.
!  totalSendMsg = 0
!  do b = 1, gr_sbNumBodies
!     bodyInfo => gr_sbBodyInfo(b)
!     if (bodyInfo % myPE == bodyInfo % bodyMaster) then
!        if (associated(bodyInfo % particlesPerProc)) then
!           do j = 1, size(bodyInfo % particlesPerProc,2)
!              if (bodyInfo % particlesPerProc(1,j) /= bodyInfo % bodyMaster) then
!                 totalSendMsg = totalSendMsg + 1
!              end if
!           end do
!           deallocate(bodyInfo % particlesPerProc) !No longer needed.
!        end if
!     end if
!  end do
!  allocate(sstatus(MPI_STATUS_SIZE,totalSendMsg))
!  allocate(sreq(totalSendMsg))
!  sreq = MPI_REQUEST_NULL

  !Count the total number of messages that we will receive so we can
  !allocate the receive status and request objects appropriately.
!  totalRecvMsg = count(gr_sbParticleCount(1:gr_sbNumBodies) > 0)
!  allocate(rstatus(MPI_STATUS_SIZE,totalRecvMsg))
!  allocate(rreq(totalRecvMsg))
!  rreq = MPI_REQUEST_NULL


  !The master process sorts the particles array and then sends
  !all off-processor particles.
!  do b = 1, gr_sbNumBodies
!     bodyInfo => gr_sbBodyInfo(b)
!     if (bodyInfo % myPE == bodyInfo % bodyMaster) then
!        call ut_sortOnProcs(totalPart, NPART_PROPS, PROC_PART_PROP, &
!             gr_meshNumProcs, bodyInfo % particles, TempBuf, &
!             perProc, logProc, countProc)
!        if(countProc>0) then
!           j=0
!           l=1
!           do i=1,countProc
!              do while(perProc(j)==0)
!                 j=j+1
!              end do
!              sendCount = perProc(j)
!              sendBufCount = sendCount*NPART_PROPS
!              if (j /= bodyInfo % bodyMaster) then
!                 if (gr_sbDebug) then
!                    write(logUnit,*) "Body", b, "master sending", &
!                         sendBufCount, " to ", j
!                 end if

!                 iSendReq = iSendReq + 1
!                 if (iSendReq > totalSendMsg) then
!                    call Driver_abortFlash("Send counter overflow")
!                 end if
!                 call MPI_Isend(bodyInfo % particles(1,l),sendBufCount,&
!                      FLASH_REAL,j,b,gr_meshComm,sreq(iSendReq),ierr)
!                 if (ierr /= MPI_SUCCESS) then
!                    call Driver_abortFlash("MPI_Isend error")
!                 endif
!              end if
!              j = j+1
!              l = l+sendCount
!           end do
!        end if
!     else
!        gettingFrom = gr_sbParticleCount(b)
!        if (gettingFrom > 0) then
!           allocate(bodyInfo % particles(NPART_PROPS,gettingFrom))
!           recvCount = gettingFrom
!           recvBufCount = recvCount * NPART_PROPS
!           if (gr_sbDebug) then
!              write(logUnit,*) "Body", b, "slave receiving", recvBufCount, &
!                   " from ", gr_sbBodyInfo(b) % bodyMaster
!           end if

!           iRecvReq = iRecvReq + 1
!           if (iRecvReq > totalRecvMsg) then
!              call Driver_abortFlash("Receive counter overflow")
!           end if
!           call MPI_Irecv(bodyInfo % particles(1,1),recvBufCount,&
!                FLASH_REAL,bodyInfo % bodyMaster,b,gr_meshComm,&
!                rreq(iRecvReq),ierr)
!           if (ierr /= MPI_SUCCESS) then
!              call Driver_abortFlash("MPI_Irecv error")
!           endif
!        end if
!     end if
!  end do

!  if (iSendReq > 0) then
!     call MPI_Waitall(iSendReq, sreq, sstatus, ierr)
!     if (ierr /= MPI_SUCCESS) then
!        call Driver_abortFlash("Send MPI_Waitall error")
!     endif
!  end if
!  if (iRecvReq > 0) then
!     call MPI_Waitall(iRecvReq, rreq, rstatus, ierr)
!     if (ierr /= MPI_SUCCESS) then
!        call Driver_abortFlash("Recv MPI_Waitall error")
!     endif
!  end if

!  call Timers_stop("particle_transfer")
!end subroutine gr_sbSendParticles


!subroutine gr_sbUpdateForces
!  use Grid_data, ONLY : gr_meshMe, gr_meshComm, gr_meshNumProcs
!  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, totalPart, &
!       gr_sbDebug, gr_sbParticleCount, solid_body
!  use Timers_interface, ONLY : Timers_start, Timers_stop
!  use ut_sortInterface, ONLY : ut_sortOnProcs
!  use Logfile_interface, ONLY : Logfile_open, Logfile_close
!  implicit none
!  type(solid_body), pointer :: bodyInfo
!  real, dimension(MDIM) :: particleposn
!  integer :: i, j, p, k, gettingFrom, ierr, blkID, b, localParticleCount, &
!       l, m, req, recvCount, logUnit, procID

!  call Timers_start("update_forces")
!  do b = 1, gr_sbNumBodies
!     bodyInfo => gr_sbBodyInfo(b)
!     if (bodyInfo % myPE == bodyInfo % bodyMaster) then
!        localParticleCount = 0
!        do i = 1, totalPart
!           if (int(bodyInfo % particles(PROC_PART_PROP,i)) == bodyInfo % bodyMaster) then
!              localParticleCount = localParticleCount + 1
!              particleposn(IAXIS) = (bodyInfo % particles(POSX_PART_PROP,i))
!              particleposn(JAXIS) = (bodyInfo % particles(POSY_PART_PROP,i))
!              particleposn(KAXIS) = (bodyInfo % particles(POSZ_PART_PROP,i))
!              blkID = int(bodyInfo % particles(BLK_PART_PROP,i))
!              call Grid_updateSolidBodyForces(blkID, b, localParticleCount, particleposn)
!           end if
!        end do
!     else
!        gettingFrom = gr_sbParticleCount(b)
!        if (gettingFrom > 0) then
!           recvCount = gettingFrom
!           do p = 1, recvCount
!              blkID = int(bodyInfo % particles(BLK_PART_PROP,p))
!              print *, "blk", RecvBuf(BLK_PART_PROP,p)
!              particleposn(IAXIS) = bodyInfo % particles(POSX_PART_PROP,p)
!              particleposn(JAXIS) = bodyInfo % particles(POSY_PART_PROP,p)
!              particleposn(KAXIS) = bodyInfo % particles(POSZ_PART_PROP,p)
!             call Grid_updateSolidBodyForces(blkID, b, p, particleposn)
              !write(*,'(a, i6, a,3f8.2,a,i6)') "myID", gr_meshMe, "receiving position", &
              !RecvBuf(POSX_PART_PROP:POSZ_PART_PROP,p), "from", gr_sbBodyInfo(b) % bodyMaster
!           enddo
!           deallocate(bodyInfo % particles)
!        end if
!     end if
!  end do
!  deallocate(gr_sbParticleCount)
!  call Timers_stop("update_forces")
!end subroutine gr_sbUpdateForces
