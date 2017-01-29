!!****if* source/physics/sourceTerms/Cool/CoolMain/Isothermal/Cool
!!
!! NAME
!!
!!  Cool(integer(IN)::blockCount
!!       integer(IN)::blockList(blockCount),
!!         real(IN) :: dt,
!!         real(IN) :: time)
!!
!!
!!
!! DESCRIPTION
!!  Apply the isothermal cooling opperator 
!!  on the list of blocks provided as input
!!
!! ARGUMENTS
!!  blockList(:) : The list of blocks on which to apply the cooling operator
!!  blockCount : The number of blocks in the list
!!  dt : the current timestep
!!  time : the current time
!!
!!***

!!REORDER(4) solnData

subroutine Cool(blockCount,blockList,dt,time)
  
  !=======================================================================
  use Cool_data, ONLY : cl_smallp, cl_smalle, cl_smallt, cl_isotherm, &
       useCool
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_getBlkIndexLimits, &
       Grid_releaseBlkPtr
  use Eos_interface, ONLY : Eos_wrapped
#include "constants.h"
#include "Flash.h"
#include "Eos.h"
  implicit none
  
  integer, intent(IN) :: blockCount
  integer,dimension(blockCount), intent(IN) :: blockList
  real, intent(IN) :: dt, time
  real, pointer, dimension(:,:,:,:) :: solnData
  integer       :: blockID, lb
  
  integer, dimension(LOW:HIGH,MDIM) :: blkLimitsGC,blkLimits
  !===========================================================================

  if(.not.useCool) return 

  do lb = 1, blockCount
     blockID=blockList(lb)
     call Grid_getBlkPtr(blockID,solnData)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     solnData(TEMP_VAR,:,:,:) = cl_isotherm
     call Eos_wrapped(MODE_DENS_TEMP, blkLimits, blockID)
     call Grid_releaseBlkPtr(blockID,solnData)
  enddo
  
  return
end subroutine Cool
