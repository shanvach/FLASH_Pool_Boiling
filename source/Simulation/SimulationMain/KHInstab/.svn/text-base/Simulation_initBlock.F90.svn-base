subroutine Simulation_initBlock(blockID)

  use Simulation_data
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getCellCoords, Grid_getDeltas, Grid_getBlkPtr, &
       Grid_releaseBlkPtr, Grid_getGeometry, Grid_renormAbundance

#include "constants.h"
#include "Flash.h"

  implicit none

  integer, intent(IN) :: blockID

  real, pointer, dimension(:,:,:,:) :: solnData
  real, allocatable, dimension(:) :: xCenter, xLeft, xRight
  real, allocatable, dimension(:) :: yCenter, yLeft, yRight
  real, allocatable, dimension(:) :: zCenter, zLeft, zRight

  integer,dimension(LOW:HIGH,MDIM)::blkLimits,blkLimitsGC
  integer :: iSize, jSize, kSize
  integer :: iSizeGC, jSizeGC, kSizeGC
  integer :: ilo, ihi
  integer :: i,j,k, istat
  real :: factor
  real :: zone, top, bot

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  iSizeGC = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  jSizeGC = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  kSizeGC = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  iSize = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
  jSize = blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
  kSize = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1

  ilo = blkLimits(LOW,IAXIS)
  ihi = blkLimits(HIGH,IAXIS)

  !! allocate all needed space
  allocate(xCenter(iSizeGC),STAT=istat)
  allocate(xLeft(iSizeGC),STAT=istat)
  allocate(xRight(iSizeGC),STAT=istat)
  allocate(yCenter(jSizeGC),STAT=istat)
  allocate(yLeft(jSizeGC),STAT=istat)
  allocate(yRight(jSizeGC),STAT=istat)
  allocate(zCenter(kSizeGC),STAT=istat)
  allocate(zLeft(kSizeGC),STAT=istat)
  allocate(zRight(kSizeGC),STAT=istat)

  xCenter(:) = 0.e0
  yCenter(:) = 0.e0
  zCenter(:) = 0.e0

  call Grid_getBlkPtr(blockID,solnData,CENTER)

  call Grid_getCellCoords(IAXIS,blockID, CENTER, .true.,xCenter,iSizeGC)
  call Grid_getCellCoords(JAXIS,blockID, CENTER, .true.,yCenter,jSizeGC)

  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              
           zone = i - blkLimits(LOW,IAXIS) + 1
           factor = 1.0 - 0.05*sin(xCenter(i)/0.05)

           top = 2.5 + 0.05*sin(PI*xCenter(i)/0.2)
           bot = 1.5 - 0.05*sin(PI*xCenter(i)/0.2)

!           if (abs(yCenter(j)) < sim_perturbRadius*factor) then
           if (yCenter(j) < top .AND. yCenter(j) > bot) then
              solnData(DENS_VAR,i,j,k) = sim_perturbDens
              solnData(VELX_VAR,i,j,k) = sim_perturbVelx
              solnData(TEMP_VAR,i,j,k) = sim_perturbTemp
              solnData(ONE_SPEC,i,j,k) = 1.0
              solnData(TWO_SPEC,i,j,k) = 0.0
           else
              solnData(DENS_VAR,i,j,k) = sim_ambientDens
              solnData(VELX_VAR,i,j,k) = sim_ambientVelx
              solnData(TEMP_VAR,i,j,k) = sim_ambientTemp
              solnData(ONE_SPEC,i,j,k) = 0.0
              solnData(TWO_SPEC,i,j,k) = 1.0
           endif

        enddo
     enddo
  enddo

  call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  
  deallocate(xLeft)
  deallocate(xRight)
  deallocate(xCenter)
  deallocate(yLeft)
  deallocate(yRight)
  deallocate(yCenter)
  deallocate(zLeft)
  deallocate(zRight)
  deallocate(zCenter)

end subroutine Simulation_initBlock
