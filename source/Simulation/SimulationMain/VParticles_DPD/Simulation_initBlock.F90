!!****if* source/Simulation/SimulationMain/VParticles_DPD/Simulation_initBlock
!!
!! NAME
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!  Simulation_initBlock(  integer(in) :: blockID,
!!                         
!!
!! DESCRIPTION   
!!     Initialize fluid properties (density, pressure, velocity, etc.) 
!!          in a single block for the unitTest/Particles
!!     Set up uniform properties with constant density and pressure
!!     For velocities, y-vel and z-vel are zero.  In the x-direction,
!!       a constant value is initialized.  Half of the y domain can have a
!!       constant multiple of the x-velocity in the rest of the domain
!!
!! ARGUMENTS
!!      blockID:     integer(in)      the current block number to be filled
!!      
!!
!! PARAMETERS
!!
!!   sim_rho_amb    Gas Density:  Entire domain receives this ambient parameter
!!   sim_p_amb      Gas Pressure: Entire domain receives this ambient parameter
!!   sim_vx_amb     Gas x-velocity:  Dominant flow velocity throughout domain 
!!   sim_vx_multiplier   Half of the domain in y has x-velocity multiplied by this value
!!
!!***


subroutine Simulation_initBlock(blockId)

!============================================================================

 
  use Simulation_data, ONLY: sim_vx_amb, sim_rho_amb, sim_p_amb, sim_vx_multiplier
  use Logfile_interface, ONLY : Logfile_stamp
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_putPointData
  implicit none 

#include "constants.h"
#include "Flash.h"


 integer, intent(IN) :: blockID

  integer               :: b, i, j, k, n, status 
  integer    :: ii, jj, kk
  integer    :: jhalf
  real       :: vx  ! velocity in x direction, may be increased in half of region with sim_vx_multiplier


  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) :: axis


  ! Write a message to stdout describing the problem setup.
  call Logfile_stamp(                              & 
          "initializing for unitTest for Particles module", '[Simulation_initBlock]')

  
!!$!----------------------------------------------------------------------------
!!$
!!$  ! get cell coordinates for the block
!!$
!!$! get the number of cells in each direction
!!$!  blkLimitsGC includes guard cells, blkLimits ignores them.
!!$  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
!!$
!!$!----------------------------------------------------------------------------
!!$
!!$! Initialize the flowfield.
!!$!  velocity in x direction is vx_amb, velocity in other directions is zero
!!$ 
!!$  ! Determine half of the y domain for multiplier
!!$  jhalf = (blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS))/2 + blkLimitsGC(LOW,JAXIS)
!!$
!!$  ! Loop over (interior) cells in the block and initialize the variables.
!!$  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
!!$     axis(KAXIS) = k
!!$
!!$     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
!!$        axis(JAXIS) = j
!!$
!!$        ! new option for changing the initial velocity with variable y
!!$        vx = sim_vx_amb
!!$        if (j .ge. jhalf) vx = sim_vx_amb*sim_vx_multiplier
!!$        
!!$
!!$        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
!!$           axis(IAXIS) = i
!!$
!!$           !-------------------------------------------------------------------
!!$           ! finally, fill the solution array
!!$           !-------------------------------------------------------------------
!!$
!!$           call Grid_putPointData(blockID, CENTER, DENS_VAR, EXTERIOR, axis, sim_rho_amb)
!!$           call Grid_putPointData(blockID, CENTER, PRES_VAR, EXTERIOR, axis, sim_p_amb)
!!$           call Grid_putPointData(blockID, CENTER, VELX_VAR, EXTERIOR, axis, vx)
!!$           call Grid_putPointData(blockID, CENTER, VELY_VAR, EXTERIOR, axis, 0.0)
!!$           call Grid_putPointData(blockID, CENTER, VELZ_VAR, EXTERIOR, axis, 0.0)
!!$
!!$
!!$        enddo  ! sweep over x
!!$     enddo     ! sweep over y
!!$  enddo        ! sweep over z
!!$
!!$  !============================================================================


  return
end subroutine Simulation_initBlock
