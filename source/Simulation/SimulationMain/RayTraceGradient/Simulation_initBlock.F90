!!****if* source/Simulation/SimulationMain/RayTraceGradient/Simulation_initBlock
!!
!! NAME
!! 
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockID)
!!                       
!!
!! DESCRIPTION
!!
!! Initial conditions for Ray Trace test problem
!!
!! ARGUMENTS
!!
!!  blockID - my block number
!!  
!!
!!***
subroutine Simulation_initBlock(blockId)

  use Simulation_data
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getCellCoords, Grid_getDeltas, Grid_getBlkPtr, &
       Grid_releaseBlkPtr, Grid_getGeometry, Grid_renormAbundance

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(IN) :: blockID
  

  real, pointer, dimension(:,:,:,:) :: solnData
  real, allocatable, dimension(:) :: xCenter, xLeft, xRight
  real, allocatable, dimension(:) :: yCenter, yLeft, yRight
  real, allocatable, dimension(:) :: zCenter, zLeft, zRight
  real, dimension(MDIM) :: delta
  real :: dx, dy, dz
  integer :: i, j, k, ii, jj, kk

  integer,dimension(LOW:HIGH,MDIM)::blkLimits,blkLimitsGC

  integer :: iSize, jSize, kSize
  integer :: iSizeGC, jSizeGC, kSizeGC
  integer :: ilo, ihi

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  iSizeGC = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  jSizeGC = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  kSizeGC = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  !! allocate all needed space
  allocate(xCenter(iSizeGC))
  allocate(yCenter(jSizeGC))
  allocate(zCenter(kSizeGC))

  xCenter(:) = 0.e0
  yCenter(:) = 0.e0
  zCenter(:) = 0.e0

  call Grid_getDeltas(blockId, delta)
  dx = delta(IAXIS)
  dy = delta(JAXIS)
  dz = delta(KAXIS)

  call Grid_getBlkPtr(blockID,solnData,CENTER)
  call Grid_getCellCoords(IAXIS,blockID, CENTER, .true.,xCenter,iSizeGC)
  call Grid_getCellCoords(JAXIS,blockID, CENTER, .true.,yCenter,jSizeGC)
  call Grid_getCellCoords(KAXIS,blockID, CENTER, .true.,zCenter,kSizeGC)

  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           
           solnData(DENS_VAR,i,j,k) = sim_gradient*sim_ncrit*(xCenter(i) + yCenter(j)) / 6.022e23
           solnData(TELE_VAR,i,j,k) = sim_temp
           solnData(EINT_VAR,i,j,k) = 0.0
           solnData(EELE_VAR,i,j,k) = 0.0

        enddo
     enddo
  enddo

  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  deallocate(xCenter)
  deallocate(yCenter)
  deallocate(zCenter)

  return
end subroutine Simulation_initBlock
