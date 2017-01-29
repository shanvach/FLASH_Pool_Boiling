!!****if* source/Particles/ParticlesMain/active/Sink/Particles_sinkComputeDt
!!
!! NAME
!!
!!  Particles_sinkComputeDt
!!
!! SYNOPSIS
!!
!!  call Particles_sinkComputeDt(integer, INTENT(in)  :: blockid,
!!                               real, INTENT(inout)  :: dt_sink,
!!                               integer, INTENT(inout)  :: dt_minloc)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   blockid : ID of block in current processor
!!
!!   dt_sink : 
!!
!!   dt_minloc : 
!!
!! NOTES
!!
!!   written by Robi Banerjee, 2007-2008
!!   modified by Christoph Federrath, 2008-2012
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!
!!***

subroutine Particles_sinkComputeDt(blockID,dt_sink,dt_minloc)

   ! Timestep correction for sink particles, idential to that
   ! of other particles. Sink particles can not travel more than
   ! a fraction sink_dt_factor across a cell in a single
   ! timestep.

   use Particles_sinkData, only : localnp, particles_local
   use Particles_data, only: pt_small, pt_globalMe
   use Grid_interface, ONLY : Grid_getBlkPhysicalSize, & 
        Grid_getBlkBoundBox, Grid_getDeltas
   use RuntimeParameters_interface, ONLY : RuntimeParameters_get

   implicit none

#include "Flash.h"
#include "constants.h"

   integer, INTENT(in)    :: blockID
   real, INTENT(inout)    :: dt_sink
   integer, INTENT(inout) :: dt_minloc(5)
   logical, save          :: first_call = .true.
   real                   :: dtx, dty, dtz, dtnew, velxabs, velyabs, velzabs
   real                   :: bsize(MDIM), delta(MDIM), lowerBound(MDIM), boundBox(2,MDIM)
   real, save             :: sink_dt_factor
   integer                :: i, lb

   if (first_call) then
      call RuntimeParameters_get("sink_dt_factor", sink_dt_factor)
      first_call =.false.
   end if

   call Grid_getBlkPhysicalSize(blockID,bsize)   ! physical size of the block in each direction
   call Grid_getBlkBoundBox(blockID,boundBox)    ! physical bounding box of the block
   lowerBound = boundBox(1,:)

   call Grid_getDeltas(blockID,delta)

   do i = 1,localnp

      lb = int(particles_local(BLK_PART_PROP,i))

      if (lb .eq. blockID) then

         dtx = HUGE(1.0)
         dty = HUGE(1.0)
         dtz = HUGE(1.0)

         velxabs = abs(particles_local(VELX_PART_PROP,i))
         if (velxabs > pt_small) then
            dtx = delta(1) / velxabs
         end if

         velyabs = abs(particles_local(VELY_PART_PROP,i))
         if (velyabs > pt_small) then
            dty = delta(2) / velyabs
         end if

         velzabs = abs(particles_local(VELZ_PART_PROP,i))
         if (velzabs > pt_small) then
            dtz = delta(3) / velzabs
         end if

         dtnew = sink_dt_factor * min(dtx, dty, dtz)
         if (dtnew .lt. dt_sink) then
            dt_sink = dtnew
            ! info about where tstep restriction took place
            dt_minloc(1) = int((particles_local(POSX_PART_PROP,i)-lowerBound(1))/delta(1)) + NGUARD+1
            dt_minloc(2) = int((particles_local(POSY_PART_PROP,i)-lowerBound(2))/delta(2)) + NGUARD+1
            dt_minloc(3) = int((particles_local(POSZ_PART_PROP,i)-lowerBound(3))/delta(3)) + NGUARD+1
            dt_minloc(4) = lb
            dt_minloc(5) = pt_globalMe
         end if

      end if     ! particle in block?

   end do   ! loop over local particles

   return

end subroutine Particles_sinkComputeDt
