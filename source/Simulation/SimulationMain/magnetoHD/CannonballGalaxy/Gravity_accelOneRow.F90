!!****if* source/Simulation/SimulationMain/magnetoHD/CannonballGalaxy/Gravity_accelOneRow
!!
!! NAME
!!
!!  Gravity_accelOneRow
!!
!!
!! SYNOPSIS
!!
!!  Gravity_accelOneRow(integer(2),intent(in):: pos, 
!!                      integer, intent(in) :: sweepDir, 
!!                      integer, intent(in) :: blockID, 
!!                      integer, intent(in) :: numCells, 
!!                      real(numCells),intent(out) :: grav, 
!!                      integer, intent(in),optional :: potentialIndex)
!!                      
!!                      
!!
!! DESCRIPTION
!!
!!  Compute components of the zone-averaged gravitational
!!  acceleration on all mesh blocks.  Either a single component
!!  of the acceleration or all three can be computed.
!!
!!  This routine computes the gravitational acceleration for a row
!!  of zones in a specified direction in a given block. First-order
!!  finite-volume differencing is used everywhere.  Formulae based
!!  on long stencils (usually high-order) may produce differences
!!  at the block boundaries for siblings as hydro solver may require
!!  several valid guard cells (e.g., PPM with parabolic
!!  interpolation for force terms needs 3 valid guard cells). Not
!!  providing such valid data may result in violation of conservation. 
!!
!! ARGUMENTS
!!
!!  pos     -       Row indices transverse to the sweep direction
!!  sweepDir   -       The sweep direction:  test against sweep_x,
!!                                 sweep_y, and sweep_z
!!  blockID   -     The local identifier of the block to work on
!!  grav()   -       Array to receive result
!!  numCells -       Number of cells to update in grav array
!!  potentialIndex      -  if specified,  Variable # to take as potential.
!!                         Default is GPOT_VAR for the potential stored in the
!!                         gpot slot of unk, which should correspond to the
!!                         potential at the current timestep.
!!
!!
!!***


subroutine Gravity_accelOneRow (pos, sweepDir, blockID, numCells, grav, potentialIndex)

  use Grid_interface, ONLY : Grid_getBlkPhysicalSize, Grid_getBlkPtr, &
    Grid_releaseBlkPtr, Grid_getCellCoords

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, dimension(2), intent(in) :: pos
  integer, intent(in)               :: sweepDir, blockID,  numCells
  real, intent(inout)               :: grav(numCells)
  integer, intent(IN),optional        :: potentialIndex
  real            :: blockSize(MDIM)
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  
  integer         :: ii, iimin, iimax, lb
  real            :: gpot(numCells), delxinv
  real, parameter :: onesixth = 1.e0/6.e0
  integer         :: potVar

  !==================================================
  
  call Grid_getBlkPhysicalSize(blockID, blockSize)
  call Grid_getBlkPtr(blockID, solnVec)
  
!! IF a variable index is explicitly specified, assume that as the potential
!! otherwise use the default current potential GPOT_VAR  
  if(present(potentialIndex)) then
     potVar=potentialIndex
  else
     potVar=GPOT_VAR
  end if

  iimin   = 1
  iimax   = numCells

  !Get row of potential values and compute inverse of zone spacing  
  if (sweepDir == SWEEP_X) then                     ! x-direction
  
     delxinv = real(NXB) / blockSize(IAXIS)
     
     gpot(:) = solnVec(potVar,:,pos(1),pos(2))

  elseif (sweepDir == SWEEP_Y) then                 ! y-direction
  
     delxinv = real(NYB) / blockSize(JAXIS)
     
     gpot(:) = solnVec(potVar,pos(1),:,pos(2))

  else                                            ! z-direction

     delxinv = real(NZB) / blockSize(KAXIS)
     
     gpot(:) = solnVec(potVar,pos(1),pos(2),:)

  endif
  
  !-------------------------------------------------------------------------------
  
  !               Compute gravitational acceleration
  
  
  !**************** first-order differences
  !                 preserves conservation
  
  delxinv = 0.5e0 * delxinv
  
  do ii = iimin+1, iimax-1
     grav(ii) = delxinv * (gpot(ii-1) - gpot(ii+1))
  enddo
  
  grav(iimin) = grav(iimin+1)     ! this is invalid data - must not be used
  grav(iimax) = grav(iimax-1)
  
  call Grid_releaseBlkPtr(blockID, solnVec)
  
  return
   
end subroutine Gravity_accelOneRow

