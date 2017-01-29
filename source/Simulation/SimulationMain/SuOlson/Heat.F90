!!****if* source/Simulation/SimulationMain/SuOlson/Heat
!!
!! NAME
!!  
!!  Heat 
!!
!!
!! SYNOPSIS
!! 
!!  call Heat (integer(IN) :: blockCount,
!!             integer(IN) :: blockList(blockCount),
!!             real(IN)    :: dt,
!!             real(IN)    :: time)
!!
!!  
!!  
!! DESCRIPTION
!!
!!  Apply the stat+gauss source term operator to a block of
!!  zones. The energy generation rate is used to update the
!!  internal energy in the zone. The phonomenological heating
!!  rate is described as a 3-D Gauss function.
!!
!!  After we call stat+gauss, call the eos to update the
!!  pressure and temperature based on the phenomenological
!!  heating.
!!  
!!
!! ARGUMENTS
!!
!!  blockCount : number of blocks to operate on
!!  blockList  : list of blocks to operate on
!!  dt         : current timestep
!!  time       : current time
!!
!!***

subroutine Heat (blockCount,blockList,dt,time)
!
!==============================================================================
!

  use Grid_interface, ONLY : Grid_getBlkPtr,         &
                             Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_releaseBlkPtr

                             
  use Eos_interface,  ONLY : Eos_wrapped
  
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  use Heat_data, ONLY : useHeat

  
#include "constants.h"
#include "Flash.h"

  implicit none
  
  integer,intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN)::blockList
  real,intent(IN) :: dt,time
  
  real    :: refTemp
  real    :: speedlt
  real    :: radconst
  real    :: ht_RadMin, ht_RadMax, radius
  real    :: ht_start, ht_end
  real    :: qdot
  real    :: gCenterX, gCenterY, gCenterZ
  integer :: i,j,k,lb,blockID
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  real,pointer,dimension(:,:,:,:):: solnData
  real,allocatable, dimension(:) :: xCenter, yCenter, zCenter
  logical :: gcell = .true.
  
  if(.not.useHeat) return  
  
  call PhysicalConstants_get("speed of light",speedlt)
  call PhysicalConstants_get("Stefan-Boltzmann",radconst)
  
  radconst = 4.0 * radconst / speedlt 
  
  !! 1000 eV
  refTemp = 1000.0*1.1604505E4
  
  ht_RadMin    = - 0.0
  ht_Radmax    = + 1.0
  
  qdot = speedlt*radconst*(refTemp)**4
  
  ht_start =  0. 
  ht_end   = 10./speedlt  
  
  gCenterX = 0.0
  gCenterY = 0.0
  gCenterZ = 0.0 
  
  if (time .ge. ht_start .and. time .le. ht_end) then
     
     do lb = 1, blockCount
        
        blockID = blockList(lb)
        
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)        
        
        allocate(xCenter (blkLimitsGC(HIGH,IAXIS)))     
        call Grid_getCellCoords(IAXIS, blockID, CENTER, gcell, & 
             xCenter, blkLimitsGC(HIGH,IAXIS))
        
#if NDIM >= 2        
        allocate(yCenter (blkLimitsGC(HIGH,JAXIS)))             
        call Grid_getCellCoords(JAXIS, blockID, CENTER, gcell, & 
             yCenter, blkLimitsGC(HIGH,JAXIS))
        
#if NDIM == 3
        allocate(zCenter (blkLimitsGC(HIGH,KAXIS)))             
        call Grid_getCellCoords(KAXIS, blockID, CENTER, gcell, & 
                   zCenter, blkLimitsGC(HIGH,KAXIS))              
#endif
#endif
        
        call Grid_getBlkPtr(blockID,solnData)        
        
        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                 
                 radius = (xCenter(i) - gCenterX)**2
                 
#if NDIM >= 2                 
                 radius = radius + (yCenter(j) - gCenterY)**2
#if NDIM == 3
                 radius = radius + (zCenter(k) - gCenterZ)**2                                  
#endif
#endif
                 radius = sqrt(radius) 
                 
                 if (radius .ge. ht_RadMin  .and. radius .le. ht_Radmax) then                                       
                    solnData(ERAD_VAR,i,j,k) = solnData(ERAD_VAR,i,j,k) + qdot*dt / solnData(DENS_VAR,  i,j,k)
                 endif
                 
              end do
           end do
        end do
        
        call Eos_wrapped(MODE_DENS_EI_GATHER,blkLimits,blockID) 
        

        deallocate(xCenter)             
#if NDIM >= 2        
        deallocate(yCenter)                     
#if NDIM == 3
        deallocate(zCenter)             
#endif
#endif     
        
     end do
     
  end if  
  
  return
  
end subroutine Heat


