!!****if* source/Simulation/SimulationMain/INavierStokes/3D/EllipseBubble_mcHYPRE/2d_ORIG_Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(in) :: blockID, 
!!                       integer(in) :: myPE)
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

subroutine Simulation_initBlock(blockID, myPE)

  use Simulation_data, ONLY : sim_xMin, sim_xMax, &
                              sim_yMin, sim_yMax, &
                              sim_gCell

  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkBoundBox,    &
                             Grid_getBlkCenterCoords

  use Driver_data, ONLY : dr_simTime

  use Multiphase_data, ONLY : mph_rho1, mph_rho2, mph_vis1, mph_vis2

  implicit none

#include "constants.h"
#include "Flash.h"

  !!$ Arguments -----------------------
  integer, intent(in) :: blockID, myPE
  !!$ ---------------------------------
 
  integer :: i, j, k
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) ::  blIndSize,blIndSizeGC

  real, dimension(MDIM)  :: coord,bsize,del
  real ::  boundBox(2,MDIM)
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData

  real :: xcell,xedge,ycell,yedge

  real :: xo, yo, rdx, rdy

  !----------------------------------------------------------------------
  
  xo = 0.5
  yo = 0.5 
  rdx = 0.35
  rdy = 0.05


  if (myPE .eq. MASTER_PE) write(*,*) 'InitBlockTime =',dr_simTime

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


  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,CENTER)
  ! Initial solution for Pressures:
  do k=1,blkLimitsGC(HIGH,KAXIS)
     do j=1,blkLimitsGC(HIGH,JAXIS)
        do i=1,blkLimitsGC(HIGH,IAXIS)

           xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
                   real(i - NGUARD - 1)*del(IAXIS) +   &
                   0.5*del(IAXIS)

           ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                   real(j - NGUARD - 1)*del(JAXIS)  +  &
                   0.5*del(JAXIS)

print*,"KPD - Simulation_initBlock... "
print*,"      DFUN is not initialized as a real "
print*,"      distance function."
print*,"      Thus the reinitialization will not "
print*,"      work, becasue the magnitude of the "
print*,"      gradient of phi cannot be = 1 !!! "

           ! Start with distance function:
           solnData(DFUN_VAR,i,j,k) = (1. - sqrt(((xcell-xo)/rdx)**2 + ((ycell-yo)/rdy)**2)) 

           ! Now densities
           if (solnData(DFUN_VAR,i,j,k) .ge. 0.0) then   ! Within ellipse bubble
              solnData(DENS_VAR,i,j,k) = mph_rho1 
              solnData(VISC_VAR,i,j,k) = mph_vis1 
           else                                          ! Outside ellipse bubble
              solnData(DENS_VAR,i,j,k) = mph_rho2
              solnData(VISC_VAR,i,j,k) = mph_vis2
           endif  

        enddo
     enddo
  enddo


  ! set values for u,v velocities and pressure
  solnData(PRES_VAR,:,:,:) = 0.0
  solnData(DELP_VAR,:,:,:) = 0.0
  solnData(DUST_VAR,:,:,:) = 0.0
  solnData(TVIS_VAR,:,:,:) = 0.0
  solnData(CURV_VAR,:,:,:) = 0.0
  solnData(SIGP_VAR,:,:,:) = 0.0  


  facexData(VELC_FACE_VAR,:,:,:) = 0.0
  faceyData(VELC_FACE_VAR,:,:,:) = 0.0
  facexData(RHDS_FACE_VAR,:,:,:) = 0.0
  faceyData(RHDS_FACE_VAR,:,:,:) = 0.0
  facexData(SIGM_FACE_VAR,:,:,:) = 0.0
  faceyData(SIGM_FACE_VAR,:,:,:) = 0.0



!!$  write(*,*) 'BlockID=',blockID
!!$  write(*,*) 'Center coordinates=',coord
!!$  write(*,*) 'Size =',bsize
!!$ 
!!$  write(*,*) 'Uveloc, 2,2,1=', &
!!$              facexData(VELC_FACE_VAR,NGUARD+1+2,NGUARD+2,1)
!!$
!!$  write(*,*) 'Vveloc, 2,2,1=', &
!!$              faceyData(VELC_FACE_VAR,NGUARD+2,NGUARD+1+2,1)


  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  call Grid_releaseBlkPtr(blockID,facexData,FACEX)
  call Grid_releaseBlkPtr(blockID,faceyData,FACEY)


  return

111    format (i4,3x,i4)
112    format (3(3x,e12.4))

end subroutine Simulation_initBlock
