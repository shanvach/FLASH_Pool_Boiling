!!****if* source/Particles/ParticlesMain/active/Sink/pt_sinkCreateParticle
!!
!! NAME
!!
!!  pt_sinkCreateParticle
!!
!! SYNOPSIS
!!
!!  call pt_sinkCreateParticle(real(in) :: x,y,z,pt
!!                             real(in) :: block_no, MyPE)
!!
!! DESCRIPTION
!!
!!   Create new particle on local CPU at position xyz
!!   accociated with block_no
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

function pt_sinkCreateParticle(x, y, z, pt, block_no, MyPE)

  use Particles_sinkData
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Particles.h"

#define get_tag(arg1,arg2) ((arg1)*65536 + (arg2))
#define get_pno(arg1) ((arg1)/65536)
#define get_ppe(arg1) ((arg1) - get_pno(arg1)*65536)

  real, intent(IN)    :: x, y, z, pt
  integer, intent(IN) :: block_no, MyPE

  logical, save       :: first_call = .true.
  integer             :: pt_sinkCreateParticle, pno, ire, tag
  logical, parameter  :: debug = .false.

  if (first_call) then
     call RuntimeParameters_get("pt_maxSinksPerProc", MaxParticlesPerProc)
     first_call = .FALSE.
  endif

  pt_sinkCreateParticle = 0

  localnp = localnp + 1
  if (localnp .GT. MaxParticlesPerProc) &
       call Driver_abortFlash('pt_sinkCreateParticle: local Particle Number exceeds MaxParticlesPerProc')

  pno = localnp

  n_empty = n_empty - 1

  do ire = 1, pt_sinkParticleProps
     particles_local(ire, pno) = 0.
  enddo

  local_tag_number = local_tag_number + 1
  
  particles_local(iptag,pno) = get_tag(local_tag_number, MyPE)
  particles_local(PROC_PART_PROP,pno) = MyPE
  particles_local(ipblk,pno) = block_no
  particles_local(ipx,pno) = x
  particles_local(ipy,pno) = y
  particles_local(ipz,pno) = z
  particles_local(ipt,pno) = pt
  particles_local(ipvx,pno) = 0.0
  particles_local(ipvy,pno) = 0.0
  particles_local(ipvz,pno) = 0.0
  particles_local(iplx,pno) = 0.0
  particles_local(iply,pno) = 0.0
  particles_local(iplz,pno) = 0.0
  particles_local(ipm,pno) = 0.0
  particles_local(ipmdot,pno) = 0.0
  particles_local(ACCX_PART_PROP,pno) = 0.0
  particles_local(ACCY_PART_PROP,pno) = 0.0
  particles_local(ACCZ_PART_PROP,pno) = 0.0
  particles_local(OACX_PART_PROP,pno) = 0.0
  particles_local(OACY_PART_PROP,pno) = 0.0
  particles_local(OACZ_PART_PROP,pno) = 0.0
  particles_local(ipdtold,pno) = 0.0

  tag = int(particles_local(iptag, pno))

  pt_sinkCreateParticle = localnp

  if (debug) then
     print*, "creating sink particle!"
     print*, "x pos=", x
     print*, "y pos=", y
     print*, "z pos=", z
     print*, "creation time=", pt
     print*, "cpu=", MyPE
     print*, "tag=",  int(particles_local(iptag, pno))
  end if

  return

end function pt_sinkCreateParticle
