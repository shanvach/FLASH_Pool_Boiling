!!****if* source/Grid/GridBoundaryConditions/gr_bcApplyToAllBlks
!!
!! NAME
!!  gr_bcApplyToAllBlks
!!
!! SYNOPSIS
!!
!!  gr_bcApplyToAllBlks(integer(IN) :: axis,
!!                      logical(IN) :: isWork)
!!  
!! DESCRIPTION 
!!
!!  This routine is a wrapper around gr_bcApplyToOneFace, and is used by UG and PM2.
!!  It calls gr_bcApplyToOneFace one each for lowerface and upperface, and repeats
!!  the process for all blocks in the grid.
!!  
!! 
!! ARGUMENTS
!!  
!!    axis           - the direction for applying BC, one of IAXIS, JAXIS, or KAXIS
!!    isWork         - is always false for UG. In PM2, if true, indicates that
!!                     the boundary conditions should be applied to work array
!!
!! NOTES
!!  A specific direction is required in axis - no ALLDIR at this point.
!!
!!***
subroutine gr_bcApplyToAllBlks_loc2D(tmpX,tmpY)
  use Grid_interface, ONLY : Grid_getLocalNumBlks, Grid_getBlkBC, &
       Grid_getBlkIndexLimits
  implicit none
#include "constants.h"
#include "Flash.h"
  
  real, dimension(GRID_IHI_GC+1,GRID_JHI_GC,GRID_KHI_GC), intent(out) :: tmpX
  real, dimension(GRID_IHI_GC,GRID_JHI_GC+1,GRID_KHI_GC), intent(out) :: tmpY

  integer,dimension(LOW:HIGH,MDIM) :: blkLimitsGC,blkLimits,blkBC
  integer :: blkNum,blockID

  blockID = 1
  call Grid_getBlkBC(blockID,blkBC)  !! For the block find the faces on the boundary
  if(blkBC(LOW,JAXIS)/=NOT_BOUNDARY) then
    tmpX(:,GRID_JLO:GRID_JLO+1,:) = 0.
    tmpY(:,GRID_JLO:GRID_JLO+1,:) = 0.
  endif
  if(blkBC(HIGH,JAXIS)/=NOT_BOUNDARY) then
    tmpX(:,GRID_JHI-1:GRID_JHI,:) = 0.
    tmpY(:,GRID_JHI-1:GRID_JHI,:) = 0.
  endif

return  
end subroutine gr_bcApplyToAllBlks_loc2D
