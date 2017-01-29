!!****h* source/physics/Flame/Flame_interface
!!
!! NAME
!!
!!  Flame_interface
!!
!! SYNOPSIS
!!
!!  use Flame_interface
!!
!! DESCRIPTION
!!
!!  Interface module for the Flame code unit.
!!
!!***

Module Flame_interface
#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  implicit none

  interface
     subroutine Flame_step( num_blocks, blockList, dt  )    
       
       integer, INTENT(in)                        :: num_blocks
       integer, INTENT(in), DIMENSION(num_blocks) :: blockList
       real,    INTENT(in)                        :: dt
       !  Evolve flame forward for all blocks in blockList by one step
       !  of size dt.
       !  Applies unsplit reaction-diffusion operater to FLAM_MSCALAR
       !  May or may not deposit energy, depending on which
       !  Flame_Effects module had been included
     end subroutine Flame_step
  end interface

  interface
     subroutine Flame_getProfile(x, f)
       
       real, intent(in)  :: x
       real, intent(out) :: f
       !  get value of progress variable a distance x from center of
       !  flame front in the steady state propagating flame (used to
       !  initialize data on mesh).
       !  x is defined such that positive x is in the direction of
       !  propagation and f = 0.5 at x = 0
     end subroutine Flame_getProfile
  end interface

  interface
     subroutine Flame_getWidth(laminarWidth)

       real, intent(OUT) :: laminarWidth
       !  approximate total width of the flame front
       !  more than this far away can initialize rpv to 0 or 1 for
       !  unburned and burned respectively
     end subroutine Flame_getWidth
  end interface
  interface
     subroutine Flame_init()

     end subroutine Flame_init
  end interface

  interface
     subroutine Flame_laminarSpeed(dens, pres, s)

       real, intent(in)   :: dens
       real, intent(in)   :: pres
       real, intent(out)  :: s
       ! return the physical laminar flame speed s
       ! behavior is strongly dependent on what FlameSpeed subunit is included
     end subroutine Flame_laminarSpeed
  end interface

  interface
     subroutine Flame_heatRelease(q, flag)

       real,                     intent(out) :: q
       integer,                  intent(in)  :: flag
       ! return how much energy the flame will release, q, in erg/gram
       ! generally for fully burning from fuel to ash
       ! flag allows for selecting burning stage for multi-stage energy
       ! release
     end subroutine Flame_heatRelease
  end interface

  interface
     subroutine Flame_rhJump(dens_u, pres_u, temp_u, ener_u, ye_u, sumy_u, &
          dens_b, pres_b, temp_b, ener_b, ye_b, sumy_b, q, s, flag)
       
       real,    intent(inout)  :: dens_u, pres_u, temp_u, ener_u
       real,    intent(out)    :: dens_b, pres_b, temp_b, ener_b
       real,    intent(in)     :: q, s
       integer, intent(in)     :: flag
       real, intent(in) ::  ye_u, sumy_u, ye_b, sumy_b
       ! Calculate the thermodynamic state of the burned material
       ! (and unburned) material by applying Rankine-Hugoniot jump condition.
       ! Unburned state is calculated frome *_u with flag as the eos mode.
       ! Burned state is calculated from energy release q (in erg/gram) and
       ! the flame speed, s, (wrt the fuel) with the composition information
       ! of the fuel specified using ye_b and sumy_b.
       ! !! This interface needs to be updated to work with or without species
     end subroutine Flame_rhJump
  end interface

  interface
     subroutine Flame_rhJumpReactive(eosData_u, qbar_u, eosData_b, qbar_b, eos_mode)

        real, dimension(EOS_NUM), intent(inout) :: eosData_u
        real,    intent(in)                     :: qbar_u
        real, dimension(EOS_NUM), intent(out)   :: eosData_b
        real,    intent(out)                    :: qbar_b
        integer, intent(in)                     :: eos_mode
        ! Calculate state of burned (and unburned) material by applying
        ! Rankine-Hugoniot jump condition.
        ! This version is for reactive ash (NSE), so that all information
        ! including the energy release (for a constant pressure burn) is
        ! calculated.
        ! Unburned state is found from calling Eos on eosData_u with mode
        ! eos_mode
        ! qbar values are in MeV/Baryon
     end subroutine Flame_rhJumpReactive
  end interface

  interface
     subroutine Flame_Effects( solnData, phi1dot, blkLimits, blkLimitsGC, time, dt, blockID)

        implicit none

        real, pointer, dimension(:,:,:,:)             :: solnData
        integer,dimension(LOW:HIGH,MDIM), intent(in)  :: blkLimits, blkLimitsGC
        real,    INTENT(in)                           :: time, dt
        integer, intent(in)                           :: blockID

#ifdef FIXEDBLOCKSIZE
        real, dimension(GRID_ILO_GC:GRID_IHI_GC,&
                        GRID_JLO_GC:GRID_JHI_GC,&
                        GRID_KLO_GC:GRID_KHI_GC), intent(in) :: phi1dot
#else
        real,dimension(:,:,:), intent(in) :: phi1dot
#endif
        ! Applies ancillary effects of the flame (e.g. energy release,
        ! composition change) by updating variables in the block passed
        ! via solnData
        ! Depending on the implementation, energy may be deposited elswhere
        ! This is really an internal fuction and could be moved to localAPI
     end subroutine Flame_Effects
  end interface

end Module Flame_interface
