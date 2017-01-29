!!****if* source/Simulation/SimulationMain/unitTest/Gravity/Poisson2/Simulation_initBlock
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
!!    Initializes fluid data (density, pressure, velocity, etc.) for
!!               a specified block.  This version sets up the Huang & Greengard
!!               Poisson solver test problem (their example 4.3).
!!
!!    Reference:   Huang, J., & Greengard, L. 2000, SIAM J. Sci. Comput.,
!!               21, 1551
!!
!!
!!
!! ARGUMENTS
!!      blockID:     integer(in)      the current block number to be filled
!!      
!!
!!
!!***



subroutine Simulation_initBlock(blockID)

  use Simulation_data, ONLY : sim_discRadius,sim_newton,sim_density
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_getBlkPtr, Grid_releaseBlkPtr

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(IN)  :: blockID
  

  real,allocatable,dimension(:) :: xCenter,yCenter,zCenter
  
  
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer       :: i, j, k, n,istat
  real :: radius, radiusSqr
  real,dimension(:,:,:,:),pointer :: solnData
  logical :: gcell=.true.
 
 
  !==============================================================


  !get the size of the block and allocate the arrays holding the coordinate info
  !this is done on a blk by block basis for compatibility with block sizes that might
  !not be fixed at compile time
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  
  sizeX = blkLimitsGC(HIGH,IAXIS)
  allocate(xCenter(sizeX),stat=istat)
  
  sizeY = blkLimitsGC(HIGH,JAXIS)
  allocate(yCenter(sizeY),stat=istat)
  
  sizeZ = blkLimitsGC(HIGH,KAXIS)
  allocate(zCenter(sizeZ),stat=istat)
  
  !initialize coord arrays
  xCenter = 0.
  yCenter = 0.
  zCenter = 0.

  
  if (NDIM > 2) then 
     call Grid_getCellCoords(KAXIS, blockID, CENTER, gcell, zCenter, sizeZ)
  endif
  
  if (NDIM > 1) then
     call Grid_getCellCoords(JAXIS, blockID, CENTER, gcell, yCenter, sizeY)
  endif
  
  call Grid_getCellCoords(IAXIS, blockID, CENTER, gcell, xCenter, sizeX)
  
  
  !  Loop over cells in the block.  For each, compute the physical
  !  position of its left and right edge and its center as well as
  !  its physical width.  Then decide whether it is inside the
  !  initial radius or outside and initialize the hydro variables
  !  appropriately.
    
  
  call Grid_getBlkPtr(blockID, solnData,CENTER)
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

           radiusSqr = xCenter(i)**2+yCenter(j)**2+zCenter(k)**2
           radius = sqrt(radiusSqr)
           if(radius < sim_discRadius)then
              solnData(DENS_VAR,i,j,k) = sim_density
              solnData(APOT_VAR,i,j,k) = 4.0*PI*sim_newton*sim_density  *  &
                   (radiusSqr/6.0 - 0.5*sim_discRadius**2)
           else
              solnData(DENS_VAR,i,j,k) = 0.0
              solnData(APOT_VAR,i,j,k) = -4.0*PI*sim_newton*sim_density *  &
                                         (sim_discRadius**3)/(3.0*radius)
           end if
        end do
     enddo
  enddo
  call Grid_releaseBlkPtr(blockID, solnData,CENTER)

  deallocate(xCenter)
  
  deallocate(yCenter)
  
  deallocate(zCenter)
  
  return

end subroutine Simulation_initBlock
