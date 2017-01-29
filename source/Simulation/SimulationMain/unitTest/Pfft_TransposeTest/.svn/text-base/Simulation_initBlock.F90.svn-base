!!****if* source/Simulation/SimulationMain/unitTest/Pfft_TransposeTest/Simulation_initBlock
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
!!  Initializes locally owned blocks to the current processor.
!! 
!! ARGUMENTS
!!
!!  blockId -          the blockId to update
!!  
!!
!!
!!***

subroutine Simulation_initBlock(blockId)
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getBlkPtr, Grid_releaseBlkPtr
  use Simulation_data, ONLY : sim_meshMe
#include "constants.h"
#include "Flash.h"

  implicit none  
  integer, intent(in) :: blockId
  

  integer :: i, j, k
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  real, pointer, dimension(:,:,:,:) :: solnData

  !I label all internal grid points with the processor ID.  
  !After all the data transfers, we use Grid_unitTest to 
  !verify that after forward and back transposes we get the same answer.
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  call Grid_getBlkPtr(blockID,solnData,CENTER)
  do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
     do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
#ifdef FLASH_GRID_UG
           solnData(DENS_VAR,i,j,k) = real(sim_meshMe)
#else
           solnData(DENS_VAR,i,j,k) = real(sim_meshMe)*100 + blockID
#endif
        end do
     end do
  end do
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)
end subroutine Simulation_initBlock
