!!****if* source/Simulation/SimulationMain/INavierStokes/2D/bhagaWeber_mcHYPRE/Simulation_initBlock
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
!!
!! 
!!
!!***

subroutine Simulation_initBlock(blockID)

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

  use Multiphase_data, ONLY : mph_rho1, mph_rho2, mph_vis1, mph_vis2, mph_meshMe

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
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData


  !----------------------------------------------------------------------
  !- kpd - Local Data
  real :: xcell,xedge,ycell,yedge
  real :: xo, yo, rdx, rdy
  real :: x, y, a, b
  real :: tsolve, t, ct, st, x1, y1, nx, ny
  integer :: iters
  real :: t1,t2,f1,f2,f,elipdist
  integer, parameter :: itermax = 50
  real, parameter :: dtmax = 1.d-8
  real, parameter :: PIHALF = 1.570796326794897d0
  !----------------------------------------------------------------------
  
  !- kpd - Bubble center coordinates 
  !xo = 0.5
  !yo = 0.5 
  xo = 0.0
  yo = 0.0 

  !- kpd - Initial bubble radii in x- and y-directions
  rdx = 0.5
  rdy = 0.5

  !if (mph_meshMe .eq. MASTER_PE) write(*,*) 'InitBlockTime =',dr_simTime

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


  !- kpd - Initialize the distance function in the 1st quadrant 
  do k=1,blkLimitsGC(HIGH,KAXIS)
     do j=1,blkLimitsGC(HIGH,JAXIS)
        do i=1,blkLimitsGC(HIGH,IAXIS)

           xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
                   real(i - NGUARD - 1)*del(IAXIS) +   &
                   0.5*del(IAXIS)

           ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                   real(j - NGUARD - 1)*del(JAXIS)  +  &
                   0.5*del(JAXIS)

          !print*,i,j,xcell,ycell

          !print*,"Cell Center",i,j,xcell,ycell

          a = rdx
          b = rdy

          !- kpd - Initialize the distance function to the circular bubble
          !solnData(DFUN_VAR,i,j,k) = (1. - sqrt(((xcell-xo)/rdx)**2 + ((ycell-yo)/rdy)**2))
          solnData(DFUN_VAR,i,j,k) = (rdx - sqrt(((xcell-xo))**2 + ((ycell-yo))**2))

         ! print*,"DFUN",i,j,k,xcell,ycell,solnData(DFUN_VAR,i,j,k)

        enddo
     enddo
  enddo

!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================


 !! Original Riaz Implementation...Start with distance function:
 !!- kpd - This was not a True distance function. 
 !solnData(DFUN_VAR,i,j,k) = (1. - sqrt(((xcell-xo)/rdx)**2 + ((ycell-yo)/rdy)**2)) 

!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================

  !- kpd - Set initial values based on distance function sign
  !----------------------------------------------------------
  do k=1,blkLimitsGC(HIGH,KAXIS)
     do j=1,blkLimitsGC(HIGH,JAXIS)
        do i=1,blkLimitsGC(HIGH,IAXIS)

           xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
                   real(i - NGUARD - 1)*del(IAXIS) +   &
                   0.5*del(IAXIS)

           ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                   real(j - NGUARD - 1)*del(JAXIS)  +  &
                   0.5*del(JAXIS)
 
           ! - kpd - Set density and viscosity initial values for each phase
           !----------------------------------------------------------------
           if (solnData(DFUN_VAR,i,j,k) .ge. 0.0) then   ! Within ellipse bubble (phase 1)
              !solnData(DENS_VAR,i,j,k) = mph_rho1 
              !solnData(VISC_VAR,i,j,k) = mph_vis1 
           else                                          ! Outside ellipse bubble (phase 2)
              !solnData(DENS_VAR,i,j,k) = mph_rho2
              !solnData(VISC_VAR,i,j,k) = mph_vis2
           endif  

        enddo
     enddo
  enddo


  !- kpd - Set initial values for velocity, pressure, surface tension, etc.
  !------------------------------------------------------------------------

  !- kpd - Cell centered variables
  solnData(PRES_VAR,:,:,:) = 0.0
  solnData(DELP_VAR,:,:,:) = 0.0
  solnData(DUST_VAR,:,:,:) = 0.0
  solnData(TVIS_VAR,:,:,:) = 0.0
  solnData(CURV_VAR,:,:,:) = 0.0
  solnData(SIGP_VAR,:,:,:) = 0.0  
  solnData(VISC_VAR,:,:,:) = 0.0 

  !- kpd - Face centered variables
  facexData(VELC_FACE_VAR,:,:,:) = 0.0
  faceyData(VELC_FACE_VAR,:,:,:) = 0.0
  facexData(RHDS_FACE_VAR,:,:,:) = 0.0
  faceyData(RHDS_FACE_VAR,:,:,:) = 0.0
  facexData(SIGM_FACE_VAR,:,:,:) = 0.0
  faceyData(SIGM_FACE_VAR,:,:,:) = 0.0


  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  call Grid_releaseBlkPtr(blockID,facexData,FACEX)
  call Grid_releaseBlkPtr(blockID,faceyData,FACEY)


  return

111    format (i4,3x,i4)
112    format (3(3x,e12.4))

end subroutine Simulation_initBlock
