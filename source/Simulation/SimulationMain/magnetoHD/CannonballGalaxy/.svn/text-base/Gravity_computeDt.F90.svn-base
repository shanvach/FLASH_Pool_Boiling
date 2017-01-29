!!****if* source/Simulation/SimulationMain/magnetoHD/CannonballGalaxy/Gravity_computeDt
!!
!! NAME
!!
!!  Gravity_computeDt
!!  
!! SYNOPSIS
!!
!!  Gravity_computeDt(integer(IN)        :: blockID,
!!                    real (OUT)         :: dt_grav,
!!                    integer(:)(INOUT)  :: dt_minloc(5))
!!
!! DESCRIPTION
!!
!!  Compute the timestep limiter due to the gravitational solver.
!!
!! ARGUMENTS
!!
!!  dt_grav:       Will Return the limiting timestep. Should be
!!                 set to a large value (1.D99) on input.
!!  dt_minloc(5):  An array to receive information about which
!!                 processor, block, and zone was responsible
!!                 for setting the limiting timestep.  The order
!!                 is i, j, k, b, p, where (i,j,k) = zone
!!                 indices, b = local block ID, and p = PE #.
!!                 This routine should only modify these values
!!                 if it changes dt_grav.
!!  blockID:       The local ID of the block to compute the
!!                 limiter on.
!!
!!***

subroutine Gravity_computeDt (blockID, dt_grav, dt_minloc)

!==============================================================================

  use Simulation_data, ONLY : sim_xCtr, sim_yCtr, sim_zCtr, sim_vrCtr, &
       sim_testSingleGalaxy, &
       sim_xMax, sim_xMin, sim_yMax, sim_yMin, sim_zMax, sim_zMin

  use Grid_interface, ONLY : Grid_getDeltas, Grid_getBlkBoundBox, &
       Grid_getMaxCommonRefinement

  use Driver_interface, ONLY : Driver_getNStep

  use Gravity_data, ONLY : grv_meshMe

  implicit none
  
#include "Flash.h"
#include "Flash_mpi.h"
#include "constants.h"

  integer, intent(IN)    ::  blockID
  integer, intent(INOUT) ::  dt_minloc(5)
  real,intent(OUT)       ::  dt_grav
  
  logical :: isOnGrid, isInBlock
  logical, save :: first_call = .true.

  real, parameter :: small = 1.0e-10
  real :: dtx, dty, dtz, dtnew
  real :: delta(MDIM), lowerBound(MDIM), boundBox(2,MDIM)

  integer, save :: lrefine_comm, nstep, nstep_last = 0
  
  integer :: lb, comm

  call Driver_getNStep(nstep)
  
  comm = MPI_COMM_WORLD

  if (first_call .or. nstep > nstep_last) then
     call Grid_getMaxCommonRefinement(comm, lrefine_comm)
     nstep_last = nstep
     first_call = .false.
  endif

  if (.not. sim_testSingleGalaxy) then

     ! Check to see if the galaxy is even on the grid first

     isOnGrid = (sim_xCtr >= sim_xMin .and. sim_xCtr <= sim_xMax)
     isOnGrid = (sim_yCtr >= sim_yMin .and. sim_yCtr <= sim_yMax) .and. isOnGrid
     isOnGrid = (sim_zCtr >= sim_zMin .and. sim_zCtr <= sim_zMax) .and. isOnGrid

     if (isOnGrid) then
          
        ! If it is, then check to see if this block has it

        call Grid_getBlkBoundBox(blockID,boundBox)    ! physical bounding box of the block
        lowerBound = boundBox(1,:)                        

        isInBlock = (sim_xCtr >= boundBox(LOW,IAXIS) .and. &
             sim_xCtr <= boundBox(HIGH,IAXIS)) .and. isInBlock
        isInBlock = (sim_yCtr >= boundBox(LOW,JAXIS) .and. &
             sim_yCtr <= boundBox(HIGH,JAXIS)) .and. isInBlock
        isInBlock = (sim_zCtr >= boundBox(LOW,KAXIS) .and. &
             sim_zCtr <= boundBox(HIGH,KAXIS)) .and. isInBlock

        if (.not. isInBlock) then
           
           ! We shouldn't be messing around in here, so bail

           return
        
        endif

        ! This is the block!

        lb = blockID

        call Grid_getDeltas(lb, delta)

     else

        ! If it's not on the grid, then we pick the size of a block at the maximum
        ! common level of refinement and say that's the size of the delta

        delta(1) = (sim_xMax - sim_xMin) / (2**(lrefine_comm-1)) / NXB
        delta(2) = (sim_yMax - sim_yMin) / (2**(lrefine_comm-1)) / NYB
        delta(3) = (sim_zMax - sim_zMin) / (2**(lrefine_comm-1)) / NZB

        ! Give a garbage value for the block

        lb = -1
        
     endif

     dtx = HUGE(1.0)
     dty = HUGE(1.0)
     dtz = HUGE(1.0)

     if (sim_vrCtr > small) dtx = delta(1) / sim_vrCtr
     if (sim_vrCtr > small) dty = delta(2) / sim_vrCtr
     if (sim_vrCtr > small) dtz = delta(3) / sim_vrCtr

     dtnew = 0.8 * min(dtx, dty, dtz)
     
     if (dtnew < dt_grav) then

        dt_grav = dtnew

        if (lb == blockID) then
           !! information about where the minimum restriction took place
           dt_minloc(1) = int((sim_xCtr-lowerBound(1))/delta(1)) + NGUARD+1
           dt_minloc(2) = int((sim_yCtr-lowerBound(2))/delta(2)) + NGUARD+1
           dt_minloc(3) = int((sim_zCtr-lowerBound(3))/delta(3)) + NGUARD+1
           dt_minloc(4) = lb
           dt_minloc(5) = grv_meshMe
        endif

     endif

  else

     ! Don't really need to worry about it

     dt_grav = huge(1.)

  endif

  
  return

end subroutine Gravity_computeDt
