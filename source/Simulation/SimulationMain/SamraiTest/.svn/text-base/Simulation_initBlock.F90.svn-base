!!****if* source/Simulation/SimulationMain/SamraiTest/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(blockId, MyPE)
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the Sod shock-tube
!!  problem.
!!
!!  Reference: Sod, G. A., 1978, J. Comp. Phys., 27, 1
!!
!!  Parameters:  blockId      The number of the block to initialize
!!
!! 
!! ARGUMENTS
!!
!!  blockId          the number of the block to update
!!
!!
!! PARAMETERS
!! 
!!  rho_left          the density to the left of the interface
!!  rho_right         the density to the right of the interface
!!
!!  p_left            the pressure to the left of the interface
!!  p_right           the pressure to the right of the interface
!!
!!  u_left            the velocity to the left of the interface
!!  u_right           the velocity to the right of the interface
!!
!!  xangle, yangle    the angle (degrees) made with respect to the x/y
!!                    axis
!!
!!  posn              the point of intersection between the shock plane
!!                    and the x-axis
!!
!!***

subroutine Simulation_initBlock(blockId)


  use Simulation_data, ONLY: sim_gCell
  use Grid_interface, ONLY : Grid_getGlobalIndexLimits, &
    Grid_getBlkIndexLimits, Grid_getBlkCornerID, Grid_putPointData

#include "constants.h"
#include "Flash.h"

  implicit none

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
  
  integer, intent(in) :: blockId
  

  integer :: i, j, k, i1, var
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) :: axis, globalRange
  integer, dimension(MDIM+1) :: dataSize 
  integer, dimension(MDIM) :: startIndex,stride
  integer :: xb,xe,yb,ye,zb,ze

  real :: tempZ,tempY,tempX,temp
  real :: pi, xPi, yPi, zPi

  integer :: beginCount
  integer, dimension(MDIM) :: position

  !in Samrai I think this is a runtime parameter?
  call Grid_getGlobalIndexLimits(globalRange)
  print *, "global Range = ", globalRange 

  pi = 8.0*atan(1.0)
  xPi = pi/globalRange(IAXIS)
  yPi = pi/globalRange(JAXIS)
  zPi = pi/globalRange(KAXIS)

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  call Grid_getBlkCornerID(blockId,startIndex,stride)
  xb = startIndex(IAXIS)
  xe = xb + blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)

  yb = startIndex(JAXIS)
  ye = yb + blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)

  zb = startIndex(KAXIS)
  ze = zb + blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)

  dataSize = 1
  axis(KAXIS) = 1
  do k = zb,ze
    ! alternate initialization values
    ! tempZ = 100*k
    tempZ = cos((k-1)*zPi)
    axis(JAXIS) = 1
    do j = yb,ye
      ! tempY = (tempZ+j)*100
      tempY = cos((j-1)*yPi)
      axis(IAXIS) = 1
      do i = xb,xe
        ! tempX = (tempY+i)*100
        ! do n = UNK_VARS_BEGIN,UNK_VARS_END
        !   temp=tempX+n
        !   call Grid_putData(axis, blockId, gcell, n, temp, dataSize)
        ! end do
        tempX = sin((i-1)*xPi)
        temp=tempX*tempY*tempZ
        do var= UNK_VARS_BEGIN,UNK_VARS_END
           call Grid_putPointData(blockId, var, beginCount, position, temp)
        end do

        axis(IAXIS) = axis(IAXIS)+1
      enddo
      axis(JAXIS) = axis(JAXIS)+1
    enddo
    axis(KAXIS) = axis(KAXIS) +1
  enddo
  
  return
end subroutine Simulation_initBlock
