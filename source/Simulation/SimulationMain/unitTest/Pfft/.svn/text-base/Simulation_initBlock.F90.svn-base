!!****if* source/Simulation/SimulationMain/unitTest/Pfft/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer (IN) ::blockId, 
!!
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes the Grid with a sinusoidal function sin(x)*cos(y)*cos(z)
!! 
!! ARGUMENTS
!!
!!  blockId -          the blockId to update
!!  
!!
!!
!!***

subroutine Simulation_initBlock(blockId)
  use Grid_interface, ONLY : Grid_getGlobalIndexLimits, &
       Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr

#include "constants.h"
#include "Flash.h"

  
  implicit none

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
  
  integer, intent(in) :: blockId
  

  integer, dimension(MDIM) :: axis, globalSize

  integer :: i, j, k, i1, var
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC

  integer, dimension(MDIM) :: blkSize
  real,pointer,dimension(:,:,:,:)::solnData
  real :: xPi, yPi, zPi, twopi

  twopi = PI*2.0

  call Grid_getGlobalIndexLimits(globalSize)

  xPi = twopi/globalSize(IAXIS)
  yPi = twopi/globalSize(JAXIS)
  zPi = twopi/globalSize(KAXIS)

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  call Grid_getBlkPtr(blockID,solnData,CENTER)
  do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           solnData(:,i,j,k)=cos(4.0*xPi*(i-1))*cos(3.0*yPi*(j-1))
        end do
     end do
  end do
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  
  return
end subroutine Simulation_initBlock
