!!****if* source/Grid/GridStructures/withTriangles/pointTopoint/bitmap/gr_sbStoreParticlesPerProc
!!
!! NAME
!!  gr_sbStoreParticlesPerProc
!!
!! SYNOPSIS
!!
!!  gr_sbStoreParticlesPerProc()
!!
!! DESCRIPTION
!!  This subroutine steps through the particles array and
!!  counts how many particles are destined for each processor.
!!  The processor ID and particle count is stored in
!!  an array in the gr_sbBodyInfo data structure.
!!
!! ARGUMENTS
!!
!!****

#include "constants.h"
#include "Flash.h"

subroutine gr_sbStoreParticlesPerProc
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_data, ONLY : gr_meshNumProcs
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, gr_sbDebug, &
       totalPart, solid_body, gr_sbParticleCount, gr_sbFirstCall
  implicit none
  integer, dimension(0:gr_meshNumProcs-1) :: perProc
  type(solid_body), pointer :: bodyInfo
  integer :: i, j, b, procID, proc, countProcs

  call Timers_start("store_particles_per_proc")
  do b = 1, gr_sbNumBodies
     bodyInfo => gr_sbBodyInfo(b)

     ! IF fixed body and not the first call:
     if ((bodyInfo%sbIsFixed .eq. CONSTANT_ONE) .and. (gr_sbFirstCall .eq. CONSTANT_ZERO)) cycle

     if (bodyInfo % myPE == bodyInfo % bodyMaster) then
        totalPart = bodyInfo % totalPart
        perProc = 0
        countProcs = 0
        do i = 1, totalPart
           procID = int(bodyInfo % particles(PROC_PART_PROP,i))
           if (procID >= 0) then
              if (perProc(procID) == 0) then
                 countProcs = countProcs + 1 !counts the no of processors that need to receive particles  
              end if
              perProc(procID) = perProc(procID) + 1
           end if
        end do

        if (countProcs > 0) then
           ! Added by Shizhao Wang on Sep 12, 2014
           if(associated(bodyInfo % particlesPerProc)) then
             nullify( bodyInfo % particlesPerProc )
           endif
           ! End of adding
           allocate(bodyInfo % particlesPerProc(1:2,1:countProcs))
           proc = 0
           do j = 0, gr_meshNumProcs-1
              if (perProc(j) > 0) then
                 proc = proc + 1
                 bodyInfo % particlesPerProc(1,proc) = j
                 bodyInfo % particlesPerProc(2,proc) = perProc(j)
              end if
           end do
        else
           nullify(bodyInfo % particlesPerProc)
        end if
     end if

     ! Set sbParticleCount to zero for this body:
     gr_sbParticleCount(b) = 0

  end do

  !We allocate and initialize in this subroutine and not 
  !gr_sbSendParticleCount to ensure the compiler flushes 
  !local initial values to memory before we start remote writes.
  !if (allocated(gr_sbParticleCount)) deallocate(gr_sbParticleCount)
  !allocate(gr_sbParticleCount(1:gr_sbNumBodies))
  !gr_sbParticleCount = 0
  call Timers_stop("store_particles_per_proc")
end subroutine gr_sbStoreParticlesPerProc
