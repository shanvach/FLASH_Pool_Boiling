!!****if* source/Grid/GridStructures/withTriangles/pointTopoint/bitmap/gr_sbUpdateForces
!!
!! NAME
!!  gr_sbUpdateForces
!!
!! SYNOPSIS
!!
!!  gr_sbUpdateForces()
!!
!! DESCRIPTION
!!
!! ARGUMENTS
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine gr_sbUpdateForces
  use Grid_data, ONLY : gr_meshMe, gr_meshComm, gr_meshNumProcs
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, totalPart, &
       gr_sbDebug, gr_sbParticleCount, solid_body
  use Timers_interface, ONLY : Timers_start, Timers_stop
  implicit none
  type(solid_body), pointer :: bodyInfo
!  real, dimension(MDIM) :: particleposn
  real :: particleData(NPART_PROPS)
  integer :: i, p, gettingFrom, blkID, b, recvCount

  call Timers_start("update_forces")
  do b = 1, gr_sbNumBodies
     bodyInfo => gr_sbBodyInfo(b)
     if (bodyInfo % myPE == bodyInfo % bodyMaster) then
        do i = 1, totalPart
           if (int(bodyInfo % particles(PROC_PART_PROP,i)) == bodyInfo % bodyMaster) then
!              particleposn(IAXIS) = (bodyInfo % particles(POSX_PART_PROP,i))
!              if(NDIM >1) then
!                 particleposn(JAXIS) = (bodyInfo % particles(POSY_PART_PROP,i))
!              endif
!              if(NDIM >2) then
!                 particleposn(KAXIS) = (bodyInfo % particles(POSZ_PART_PROP,i))
!              endif
              blkID = int(bodyInfo % particles(BLK_PART_PROP,i))
              particleData = bodyInfo % particles(1:NPART_PROPS,i)
              call Grid_updateSolidBodyForces(blkID, particleData)
           end if
        end do
     else
        gettingFrom = gr_sbParticleCount(b)
        if (gettingFrom > 0) then
           recvCount = gettingFrom
           do p = 1, recvCount
              blkID = int(bodyInfo % particles(BLK_PART_PROP,p))
!              particleposn(IAXIS) = bodyInfo % particles(POSX_PART_PROP,p)
!              if(NDIM >1) then
!                 particleposn(JAXIS) = bodyInfo % particles(POSY_PART_PROP,p)
!              endif
!              if(NDIM >2) then
!                 particleposn(KAXIS) = bodyInfo % particles(POSZ_PART_PROP,p)
!              endif
              particleData = bodyInfo % particles(1:NPART_PROPS,p)
              call Grid_updateSolidBodyForces(blkID, particleData)
              !write(*,'(a, i6, a,3f8.2,a,i6)') "myID", gr_meshMe, "receiving position", & 
              !RecvBuf(POSX_PART_PROP:POSZ_PART_PROP,p), "from", gr_sbBodyInfo(b) % bodyMaster 
           enddo
!           deallocate(bodyInfo % particles)
        end if
     end if
  end do
!  deallocate(gr_sbParticleCount)
  call Timers_stop("update_forces")
end subroutine gr_sbUpdateForces
