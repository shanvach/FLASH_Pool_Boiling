!!****if* source/Particles/ParticlesMain/active/Sink/pt_sinkDumpParticles
!!
!! NAME
!!
!!  pt_sinkDumpParticles
!!
!! SYNOPSIS
!!
!!  call pt_sinkDumpParticles(real, intent(IN)  :: simtime)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   simtime : 
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

!
! writes all important sink particle properties (positions, ...) each timestep to file
!
subroutine pt_sinkDumpParticles(simtime)

  use Particles_sinkData
  use pt_sinkInterface, only: pt_sinkGatherGlobal
  use Driver_data, ONLY : dr_globalMe

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Particles.h"

  real, intent(IN)   :: simtime

  integer, parameter :: funit_evol = 15
  character(len=80)  :: outfile = "sinks_evol.dat"
  integer            :: i
  integer, save      :: MyPE, MasterPE, ipax, ipay, ipaz
  logical, save      :: firstCall = .TRUE.

  if (firstCall) then

     MyPE     = dr_globalMe
     MasterPE = MASTER_PE
     ipax     = ACCX_PART_PROP
     ipay     = ACCY_PART_PROP
     ipaz     = ACCZ_PART_PROP

     if (MyPE .eq. MasterPE) then
        open(funit_evol, file=trim(outfile), position='APPEND')
        write(funit_evol,'(20(1X,A16))') '[00]part_tag', '[01]time', '[02]posx', '[03]posy', '[04]posz', &
                                         '[05]velx', '[06]vely', '[07]velz', '[08]accelx', '[09]accely', &
                                         '[10]accelz', '[11]anglx', '[12]angly', '[13]anglz', '[14]mass', &
                                         '[15]mdot', '[16]mgas', '[17]ptime'
        close(funit_evol)
     endif

     firstCall = .false.

  endif

  ! exchange particle information across CPUs (this needs to be called by all processors)
  call pt_sinkGatherGlobal()

  ! only the master processor dumps the data of all particles (global list perticlest)
  if (MyPE .NE. MasterPE) return

  open(funit_evol, file=trim(outfile), position='APPEND')
  do i = 1, localnpf

     write(funit_evol,'(1(1X,I16),19(1X,ES16.9))') &
          int(particles_global(iptag,i)), &
          simtime, &
          particles_global(ipx,i), &
          particles_global(ipy,i), &
          particles_global(ipz,i), &
          particles_global(ipvx,i), &
          particles_global(ipvy,i), &
          particles_global(ipvz,i), &
          particles_global(ipax,i), &
          particles_global(ipay,i), &
          particles_global(ipaz,i), &
          particles_global(iplx,i), &
          particles_global(iply,i), &
          particles_global(iplz,i), &
          particles_global(ipm,i), &
          particles_global(ipmdot,i), &
          particles_global(ipmgas,i), &
          particles_global(ipt,i)

  enddo
  close(funit_evol)

  return

end subroutine pt_sinkDumpParticles
