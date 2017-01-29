!!****if* source/Grid/GridStructures/withTriangles/pointTopoint/bitmap/gr_sbSendParticleCount
!!
!! NAME
!!  gr_sbSendParticleCount
!!
!! SYNOPSIS
!!
!!  gr_sbSendParticleCount()
!!
!! DESCRIPTION
!!
!!  Use the processor ID and particle count information
!!  collected in gr_sbStoreParticlesPerProc subroutine to inform 
!!  the required slave processors how many particles they will
!!  receive. The communication depends upon MPI-2 and is
!!  completely contained between two synchronization fences.
!!
!! ARGUMENTS
!!
!! ***

#include "constants.h"
#include "Flash.h"

subroutine gr_sbSendParticleCount
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_data, ONLY : gr_meshComm
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, gr_sbDebug, &
       solid_body, gr_sbParticleCount, gr_sbFirstCall
  implicit none
  include "Flash_mpi.h"
  type(solid_body), pointer :: bodyInfo
  integer :: win, ierr, b, j, targetRank, dispUnit
  integer (kind=MPI_ADDRESS_KIND) :: lowerbound, intSize, winSize, targetDisp

  call Timers_start("send_particle_count")
  call MPI_Type_get_extent(MPI_INTEGER, lowerbound, intSize, ierr)
  winSize = intSize * gr_sbNumBodies
  dispUnit = intSize
  call MPI_Win_create(gr_sbParticleCount, winSize, dispUnit, &
       MPI_INFO_NULL, gr_meshComm, win, ierr)

  call MPI_Win_fence(0, win, ierr)
  do b = 1, gr_sbNumBodies

     ! IF fixed body and not the first call:
     if ((gr_sbBodyInfo(b)%sbIsFixed .eq. CONSTANT_ONE) .and. (gr_sbFirstCall.eq. CONSTANT_ZERO)) cycle

     bodyInfo => gr_sbBodyInfo(b)
     if (bodyInfo % myPE == bodyInfo % bodyMaster) then
        if (associated(bodyInfo % particlesPerProc)) then
           do j = 1, size(bodyInfo % particlesPerProc,2)
              targetRank = bodyInfo % particlesPerProc(1,j)
              if (targetRank /= bodyInfo % bodyMaster) then
                 !print *, "PE", bodyInfo % bodyMaster, "sending", & 
                 !     bodyInfo % particlesPerProc(2,j), "to", &
                 !     targetRank, "for body:", b 
                 targetDisp = b-1
                 !Note that MPI_Put behaves like MPI_Isend and so 
                 !communication may be deferred until the MPI_Win_fence.                                            
                 call MPI_Put(bodyInfo % particlesPerProc(2,j), 1, MPI_INTEGER, &
                      targetRank, targetDisp, 1, MPI_INTEGER, win, ierr)
              end if
           end do
        end if
     end if

  end do
  call MPI_Win_fence(0, win, ierr)

  call MPI_Win_free(win, ierr)
  call Timers_stop("send_particle_count")
end subroutine gr_sbSendParticleCount
