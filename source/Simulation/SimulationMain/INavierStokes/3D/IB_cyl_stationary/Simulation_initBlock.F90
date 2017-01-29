!!****if* source/Simulation/SimulationMain/INavierStokes/2D/LidDrivenCavity/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(in) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.
!!
!!  Reference:
!!
!! 
!! ARGUMENTS
!!
!!  blockID -          the number of the block to update
!!  myPE   -           my processor number
!!
!! 
!!
!!***

subroutine Simulation_initBlock(blockId)

  use Simulation_data, ONLY : sim_xMin, sim_xMax, &
                              sim_yMin, sim_yMax, &
                              sim_zMin, sim_zMax, &
                              sim_gCell

  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkBoundBox,    &
                             Grid_getBlkCenterCoords

  use Driver_data, ONLY : dr_simTime

  implicit none

#include "constants.h"
#include "Flash.h"

  !!$ Arguments -----------------------
  integer, intent(in) :: blockID
  !!$ ---------------------------------
 
  integer :: i, j, k
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) ::  blIndSize,blIndSizeGC

  real, dimension(MDIM)  :: coord,bsize,del
  real ::  boundBox(2,MDIM)
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData


  real :: xedge(GRID_IHI_GC+1),xcell(GRID_IHI_GC)
  real :: yedge(GRID_JHI_GC+1),ycell(GRID_JHI_GC)
  real :: zedge(GRID_KHI_GC+1),zcell(GRID_KHI_GC)

  real, parameter :: Uoo = 1.
  real, parameter :: a = 0.5

  real :: r,theta,vr,vtheta

  !----------------------------------------------------------------------
  
  !if (myPE .eq. MASTER_PE) write(*,*) 'InitBlockTime =',dr_simTime

  ! Get nxb, nyb and nxb:
  !call Grid_getBlkIndexSize(blockId,blIndSize,blIndSizeGC)

  !nxb = blIndSize(1)
  !nyb = blIndSize(2)
  !nzb = blIndSize(3)

  ! Get Coord and Bsize for the block:
  ! Bounding box:
  call Grid_getBlkBoundBox(blockId,boundBox)
  bsize(:) = boundBox(2,:) - boundBox(1,:)

  call Grid_getBlkCenterCoords(blockId,coord)

  ! Get blocks dx, dy ,dz:
  call Grid_getDeltas(blockID,del)

  ! Point to Blocks centered variables:
  call Grid_getBlkPtr(blockID,solnData,CENTER)

  ! Point to Blocks face variables: 
  call Grid_getBlkPtr(blockID,facexData,FACEX)
  call Grid_getBlkPtr(blockID,faceyData,FACEY)
  call Grid_getBlkPtr(blockID,facezData,FACEZ)

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,CENTER)

  ! Compute Grid line locations:
  ! Z:
  do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     zedge(k) = coord(KAXIS) - bsize(KAXIS)/2.0 + real(k - NGUARD - 1)*del(KAXIS)
     zcell(k) = zedge(k) + 0.5*del(KAXIS)     
  enddo 
  zedge(blkLimitsGC(HIGH,KAXIS)+1)=coord(KAXIS)-bsize(KAXIS)/2.0 + &
                                   real(blkLimitsGC(HIGH,KAXIS)-NGUARD)*del(KAXIS)
  ! Y:
  do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
     yedge(j) = coord(JAXIS) - bsize(JAXIS)/2.0 + real(j - NGUARD - 1)*del(JAXIS)
     ycell(j) = yedge(j) + 0.5*del(JAXIS)
  enddo
  yedge(blkLimitsGC(HIGH,JAXIS)+1)=coord(JAXIS)-bsize(JAXIS)/2.0 + &
                                   real(blkLimitsGC(HIGH,JAXIS)-NGUARD)*del(JAXIS)

  ! X:
  do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
     xedge(i) = coord(IAXIS) - bsize(IAXIS)/2.0 + real(i - NGUARD - 1)*del(IAXIS)
     xcell(i) = xedge(i) + 0.5*del(IAXIS)
  enddo
  xedge(blkLimitsGC(HIGH,IAXIS)+1)=coord(IAXIS)-bsize(IAXIS)/2.0+ &
                                   real(blkLimitsGC(HIGH,IAXIS)-NGUARD)*del(IAXIS)

  facexData(VELC_FACE_VAR,:,:,:) = 0.0
  !faceyData(VELC_FACE_VAR,:,:,:) = 0.0
  faceyData(VELC_FACE_VAR,:,:,:) = 1.
  facezData(VELC_FACE_VAR,:,:,:) = 0.0

#ifdef POTENTIAL_FLOW
  ! Velocities: Cylinder center assumed in {y,z}={0,0} for all x.
  ! Z velocity:
  do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)+1
     do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)  
          theta = atan2(ycell(j),zedge(k))
          r     = sqrt(ycell(j)**2. + zedge(k)**2.)
          if ( r .gt. a) then
            vr     = Uoo*(1.-a**2./r**2.)*cos(theta)
            vtheta =-Uoo*(1.+a**2./r**2.)*sin(theta) 
            facezData(VELC_FACE_VAR,i,j,k) = vr*cos(theta) - vtheta*sin(theta)
          endif
        enddo
     enddo
  enddo

  ! Y velocity:
  do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)+1
        do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
          theta = atan2(yedge(j),zcell(k))
          r     = sqrt(yedge(j)**2. + zcell(k)**2.)
          if ( r .gt. a) then
            vr     = Uoo*(1.-a**2./r**2.)*cos(theta)
            vtheta =-Uoo*(1.+a**2./r**2.)*sin(theta)
            faceyData(VELC_FACE_VAR,i,j,k) = vr*sin(theta) + vtheta*cos(theta)
          endif
        enddo
     enddo
  enddo
#endif

  ! set values for u,v velocities and pressure
  solnData(PRES_VAR,:,:,:) = 0.0
  solnData(DELP_VAR,:,:,:) = 0.0
  solnData(DUST_VAR,:,:,:) = 0.0
  solnData(TVIS_VAR,:,:,:) = 0.0

  facexData(RHDS_FACE_VAR,:,:,:) = 0.0
  faceyData(RHDS_FACE_VAR,:,:,:) = 0.0
  facezData(RHDS_FACE_VAR,:,:,:) = 0.0

  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  call Grid_releaseBlkPtr(blockID,facexData,FACEX)
  call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
  call Grid_releaseBlkPtr(blockID,facezData,FACEZ)



  return

111    format (i4,3x,i4)
112    format (3(3x,e12.4))

end subroutine Simulation_initBlock
