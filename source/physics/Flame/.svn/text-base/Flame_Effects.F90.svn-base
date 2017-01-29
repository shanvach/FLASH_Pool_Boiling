!!****f* source/physics/Flame/Flame_Effects
!!
!! NAME
!!
!!  Flame_Effects
!!
!! SYNOPSIS
!!
!!  call Flame_Effects( real(:,:,:,:) :: solnData,
!!                      real(:,:,:) :: phi1dot,
!!                      integer(:,:) :: blkLimits,
!!                      integer(:,:) :: blkLimitsGC,
!!                      real :: time,
!!                      real :: dt,
!!                      integer :: blockID )
!!
!! DESCRIPTION
!!
!!   Uses phi1dot to update energy and local material properties based on flame progress
!!
!! ARGUMENTS
!!
!!   blockID -    ID of the block from to whose data solnData is pointing.
!!
!!***

!!REORDER(4): solnData
#ifdef DEBUG_ALL
!#define DEBUG_INTERFACES
#endif

subroutine Flame_Effects( solnData, phi1dot, blkLimits, blkLimitsGC, time, dt, blockID)

  implicit none

#include "constants.h"
#include "Flash.h"
    
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

  return
end subroutine Flame_Effects

