!!****if* source/Simulation/SimulationMain/INavierStokes/2D/IB_Cyl_parallelIBVP_HYPRE_VD_halfDiamBel/Simulation_initBlock
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
                              sim_gCell, sim_waveA, &
                              sim_nuc_site_x, sim_nuc_site_y,&
                              sim_nuc_site_z, sim_nuc_radii,&
                              sim_nucSiteDens, sim_Tbulk

  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_getCellMetrics,    &
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

  real, dimension(GRID_IHI_GC) :: delx
  real, dimension(GRID_JHI_GC) :: dely
  real, dimension(GRID_KHI_GC) :: delz

  real, dimension(MDIM)  :: coord,bsize,del
  real ::  boundBox(2,MDIM)
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData

  real :: xcell,xedge,ycell,yedge,zcell

  real, dimension(GRID_IHI_GC,3,1) :: dx
  real, dimension(GRID_JHI_GC,3,1) :: dy
  real, dimension(GRID_KHI_GC,3,1) :: dz

  real :: xCoord(GRID_IHI_GC), xFace(GRID_IHI_GC+1), yCoord(GRID_JHI_GC)
 

  real :: A0

  real :: A,B,emp,fs,x0,y0,r0,solnX,z0,x1,y1,z1,x2,y2,z2,d1,d2,d3

  real :: x3,y3,z3,d4,r1

  real :: fn(10)
  real :: x4,x5,x6,x7,x8
  real :: y4,y5,y6,y7,y8
  real :: d5,d6,d7,d8,d9
  real :: z4,z5,z6,z7,z8 

  real, dimension(17) :: Nuc_radii,Nuc_sites_x,Nuc_sites_z,Nuc_sites_y
  real :: Nuc_dfun
  integer :: Nuc_Index, bli
  real :: th_radii,centerx,centery,radius
  real :: xl,xr,yl,yr,dxl,dxr,dyl,dyr

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
  !call Grid_getDeltas(blockID,del)

  call Grid_getCellCoords(IAXIS,blockID,CENTER,.true.,xCoord,GRID_IHI_GC)
  call Grid_getCellCoords(JAXIS,blockID,CENTER,.true.,yCoord,GRID_JHI_GC)
  call Grid_getCellCoords(IAXIS,blockID,FACES,.true.,xFace,GRID_IHI_GC+1) 

  call Grid_getCellMetrics(IAXIS,blockID,LEFT_EDGE, .true.,dx(:,LEFT_EDGE,blockID), GRID_IHI_GC) 
  call Grid_getCellMetrics(IAXIS,blockID,CENTER,    .true.,dx(:,CENTER,blockID),    GRID_IHI_GC) 
  call Grid_getCellMetrics(IAXIS,blockID,RIGHT_EDGE,.true.,dx(:,RIGHT_EDGE,blockID),GRID_IHI_GC) 

  call Grid_getCellMetrics(JAXIS,blockID,LEFT_EDGE, .true.,dy(:,LEFT_EDGE,blockID), GRID_JHI_GC) 
  call Grid_getCellMetrics(JAXIS,blockID,CENTER,    .true.,dy(:,CENTER,blockID),    GRID_JHI_GC) 
  call Grid_getCellMetrics(JAXIS,blockID,RIGHT_EDGE,.true.,dy(:,RIGHT_EDGE,blockID),GRID_JHI_GC) 

#if NDIM ==3
  call Grid_getCellMetrics(KAXIS,blockID,LEFT_EDGE, .true.,dz(:,LEFT_EDGE,blockID), GRID_KHI_GC) 
  call Grid_getCellMetrics(KAXIS,blockID,CENTER,    .true.,dz(:,CENTER,blockID),    GRID_KHI_GC) 
  call Grid_getCellMetrics(KAXIS,blockID,RIGHT_EDGE,.true.,dz(:,RIGHT_EDGE,blockID),GRID_KHI_GC) 
#endif  

 ! Point to Blocks centered variables:
  call Grid_getBlkPtr(blockID,solnData,CENTER)

  ! Point to Blocks face variables: 
  call Grid_getBlkPtr(blockID,facexData,FACEX)
  call Grid_getBlkPtr(blockID,faceyData,FACEY)

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,CENTER)

  centerx = sim_xMin+(sim_xMax-sim_xMin)/2
  centery = sim_yMin+(sim_yMax-sim_yMin)/4
  radius = 0.5
  !- kpd - Initialize the distance function in the 1st quadrant 
  do k=1,blkLimitsGC(HIGH,KAXIS)
     do j=1,blkLimitsGC(HIGH,JAXIS)
        do i=1,blkLimitsGC(HIGH,IAXIS)

           xcell = xCoord(i)

           ycell = yCoord(j)

           zcell = 0.0

           solnData(DFUN_VAR,i,j,k) = -(sqrt((xcell-centerx)**2+(ycell-centery)**2+zcell**2)-radius)!,ycell+23.0)

        enddo
     enddo
  enddo


  !- kpd - Initialize the distance function in the 1st quadrant 
  do k=1,blkLimitsGC(HIGH,KAXIS)
     do j=1,blkLimitsGC(HIGH,JAXIS)+1
        do i=1,blkLimitsGC(HIGH,IAXIS)

          faceyData(VELC_FACE_VAR,i,j,k) = 0.0 !1.5*(1.0 - (xcell*xcell)/(cheight*cheight))

        enddo
     enddo
  enddo

#if(0)
  !- wsz - Initialize the velocity in the 1st quadrant
 
  do k=1,blkLimitsGC(HIGH,KAXIS) 
     do j=1,blkLimitsGC(HIGH,JAXIS)
        do i=1,blkLimitsGC(HIGH,IAXIS)+1

           ycell  = yCoord(j)

          if (ycell .LE. 0.0) then
             facexData(VELC_FACE_VAR,i,j,k) = 0.0d0
          else
             facexData(VELC_FACE_VAR,i,j,k) = 0.0d0
          end if

        enddo
     enddo
  enddo
#endif

  ! set values for u,v velocities and pressure
  solnData(PRES_VAR,:,:,:) = 0.0
  solnData(DELP_VAR,:,:,:) = 0.0
  solnData(DUST_VAR,:,:,:) = 0.0
  solnData(TVIS_VAR,:,:,:) = 0.0 

  solnData(POOD_VAR,:,:,:) = 0.0
  facexData(SIGO_FACE_VAR,:,:,:) = 0.0
  facexData(SIGOO_FACE_VAR,:,:,:) = 0.0
  faceyData(SIGO_FACE_VAR,:,:,:) = 0.0
  faceyData(SIGOO_FACE_VAR,:,:,:) = 0.0

  solnData(CURV_VAR,:,:,:) = 0.0
  solnData(SIGP_VAR,:,:,:) = 0.0
  solnData(VISC_VAR,:,:,:) = 0.0
  solnData(PFUN_VAR,:,:,:) = 0.0

  facexData(RHDS_FACE_VAR,:,:,:) = 0.0
  faceyData(RHDS_FACE_VAR,:,:,:) = 0.0
  facexData(SIGC_FACE_VAR,:,:,:) = 0.0
  faceyData(SIGC_FACE_VAR,:,:,:) = 0.0

  facexData(SIGM_FACE_VAR,:,:,:) = 0.0
  faceyData(SIGM_FACE_VAR,:,:,:) = 0.0
  facexData(RH1F_FACE_VAR,:,:,:) = 0.0
  faceyData(RH1F_FACE_VAR,:,:,:) = 0.0
  facexData(RH2F_FACE_VAR,:,:,:) = 0.0
  faceyData(RH2F_FACE_VAR,:,:,:) = 0.0

  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  call Grid_releaseBlkPtr(blockID,facexData,FACEX)
  call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
  return

111    format (i4,3x,i4)
112    format (3(3x,e12.4))

end subroutine Simulation_initBlock
