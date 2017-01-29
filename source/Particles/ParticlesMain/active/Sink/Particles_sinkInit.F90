!!****if* source/Particles/ParticlesMain/active/Sink/Particles_sinkInit
!!
!! NAME
!!
!!  Particles_sinkInit
!!
!! SYNOPSIS
!!
!!  call Particles_sinkInit(logical, INTENT(in)  :: restart)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   restart : 
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

subroutine Particles_sinkInit ( restart)

   ! initialize the sink particles arrays

   use Particles_data, ONLY : pt_globalMe
   use Particles_sinkData
   use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
        RuntimeParameters_mapStrToInt
   use Driver_interface, ONLY : Driver_abortFlash

   implicit none

#include "constants.h"
#include "Flash.h"
#include "Particles.h"

   logical, INTENT(in) :: restart

   integer :: ierr

   call RuntimeParameters_get ("pt_maxSinksPerProc", pt_maxSinksPerProc)
   call RuntimeParameters_get ("UseSinkParticles", useSinkParticles)

   if (.not. UseSinkParticles) return

   allocate (NumParticlesPerBlock(MAXBLOCKS))
   NumParticlesPerBlock(:) = 0

   if (pt_globalMe .eq. MASTER_PE) print*, "Initializing Sink Particles"

   ! Sink particle Properties (see Flash.h)
   ipx = POSX_PART_PROP
   ipy = POSY_PART_PROP
   ipz = POSZ_PART_PROP
   ipvx = VELX_PART_PROP
   ipvy = VELY_PART_PROP
   ipvz = VELZ_PART_PROP
   ipm = MASS_PART_PROP
   ipblk = BLK_PART_PROP

   iptag = TAG_PART_PROP
   iplx = X_ANG_PART_PROP
   iply = Y_ANG_PART_PROP
   iplz = Z_ANG_PART_PROP
   ipt = CREATION_TIME_PART_PROP
   ipmdot = ACCR_RATE_PART_PROP 
   ipraccr = ACCR_RADIUS_PART_PROP
   iold_pmass = OLD_PMASS_PART_PROP
   ipmgas = MGAS_PART_PROP
   ipdtold = DTOLD_PART_PROP
   ipbflx = X_BFLUX_PART_PROP
   ipbfly = Y_BFLUX_PART_PROP
   ipbflz = Z_BFLUX_PART_PROP

   n_empty = MaxParticlesPerProc
   RunningParticles = .true.
   if (.not. restart) then
      localnp = 0
      localnpf = 0
   end if

   if (.not. restart) then !if we starting from scratch

      if (.not. allocated(particles_local)) then
         allocate (particles_local(pt_sinkParticleProps,pt_maxSinksPerProc), stat=ierr)
         if (ierr /= 0) then
            call Driver_abortFlash("Particles_init:  could not allocate particles_local array")
         endif
      endif

      if (.not. allocated(particles_global)) then
           allocate (particles_global(pt_sinkParticleProps,pt_maxSinksPerProc), stat=ierr)
           if (ierr /= 0) then
              call Driver_abortFlash("Particles_init:  could not allocate particles_global array for sink particles")
           endif
        end if

      particles_local = NONEXISTENT
      particles_global = NONEXISTENT

   end if  ! end of .not. restart

   if (allocated(particles_local)) particles_local(VELX_PART_PROP,:)=0.0
   if (allocated(particles_global)) particles_global(VELX_PART_PROP,:)=0.0

   if (allocated(particles_local)) particles_local(VELY_PART_PROP,:)=0.0
   if (allocated(particles_global)) particles_global(VELY_PART_PROP,:)=0.0

   if (allocated(particles_local)) particles_local(VELZ_PART_PROP,:)=0.0
   if (allocated(particles_global)) particles_global(VELZ_PART_PROP,:)=0.0

   return

end subroutine Particles_sinkInit
