!!****if* source/Simulation/SimulationMain/unitTest/Particles/sim_randomMove
!!
!! NAME
!!
!!  sim_randomMove
!!
!! SYNOPSIS
!!
!!  sim_randomMove( integer(in) :: blockCount,
!!                         integer(in) :: blockList(:) )
!!
!! ARGUMENTS
!!
!!   blockCount:     IN    integer  Number of blocks on the current processor 
!!   blockList(:):   IN    integer  list of LEAF blocks 
!!  
!! DESCRIPTION
!!  
!!  Apply a random increment to the velocity field.  This routine replaces
!!    the standard Hydro/physics in order to isolate particle movement.
!!
!! PARAMETERS
!!
!!   sim_vx_amb     Gas x-velocity:  Dominant flow velocity throughout domain 
!!   sim_vx_multiplier   Half of the domain in y has x-velocity multiplied by this value
!!   sim_seed   Random number seed -- NOT USED please ignore
!!   sim_vx_pert   Scales [-1,1] random number in x direction: set to zero for uniform flow
!!   sim_vy_pert   Scales [-1,1] random number in y direction: set to zero for uniform flow
!!   sim_vz_pert   Scales [-1,1] random number in z direction: set to zero for uniform flow
!!
!!***

subroutine sim_randomMove(blockCount,blockList)
!===============================================================================

  use Simulation_data, ONLY: sim_vx_pert, sim_vy_pert, sim_vz_pert, sim_seed
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getPointData, Grid_putPointData
  use Particles_data, ONLY:
  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(IN) :: blockCount
  integer, intent(IN) :: blockList(blockCount)

  real           ::  velX, velY, velZ
  real           ::  randomX, randomY, randomZ

  integer        :: b, i, j, k, seed_size
  integer        :: blockID
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM)   :: axis
  integer, dimension(10)     :: seed_array


!===============================================================================
  
  ! loop over all blocks
  do b = 1, blockCount
     blockID = blockList(b)

     !----------------------------------------------------------------------------
     ! get the number of cells in each direction for the block
     !  blkLimitsGC includes guard cells, blkLimits ignores them.
     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)

     !----------------------------------------------------------------------------

     !  seed the random number generator for reproducibility
!  D'Oh!  this doesn't produce random numbers over each time step
!     seed_array = sim_seed
!     call random_seed( put=seed_array(1:10))

     ! Loop over (interior) cells in the block and initialize the variables.


     do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
        axis(KAXIS) = k

        do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
           axis(JAXIS) = j

           do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              axis(IAXIS) = i

              ! Generate random numbers between 0 and 1
              ! and rescale to -1 to 1
              call random_number(randomX)
              randomX = randomX * 2.0 - 1.0
              if (NDIM >= 2) then
                 call random_number(randomY)
                 randomY = randomY * 2.0 - 1.0
              else
                 randomY = 0.0
              endif
              if (NDIM == 3) then
                 call random_number(randomZ)
                 randomZ = randomZ * 2.0 - 1.0
              else
                 randomZ = 0.0
              endif
!!$              if (i==1.and.j==1) print *,'random numbers',randomX,randomY,randomZ

              ! scale by the user requested perturbation
              randomX = randomX * sim_vx_pert
              randomY = randomY * sim_vy_pert
              randomZ = randomZ * sim_vz_pert
!!$              if (i==1.and.j==1) print *,'perturbed random',randomX,randomY,randomZ

              !  Get the velocity data
              call Grid_getPointData(blockID, CENTER, VELX_VAR, EXTERIOR, axis, velX)
              call Grid_getPointData(blockID, CENTER, VELY_VAR, EXTERIOR, axis, velY)
              call Grid_getPointData(blockID, CENTER, VELZ_VAR, EXTERIOR, axis, velZ)

              ! increment by the random perturbation
!!$              if (i==1.and.j==1) print *,'original velocities',velX,velY,velZ
              velX = velX + randomX
              velY = velY + randomY
              velZ = velZ + randomZ
!!$              if (i==1.and.j==1) print *,'modified velocities',velX,velY,velZ


              ! finally, fill the solution array
              call Grid_putPointData(blockID, CENTER, VELX_VAR, EXTERIOR, axis, velX)
              call Grid_putPointData(blockID, CENTER, VELY_VAR, EXTERIOR, axis, velY)
              call Grid_putPointData(blockID, CENTER, VELZ_VAR, EXTERIOR, axis, velZ)


           enddo  ! sweep over x
        enddo     ! sweep over y
     enddo        ! sweep over z
  enddo          ! loop over blocks



end subroutine sim_randomMove
