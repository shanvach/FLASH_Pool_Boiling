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
                              sim_gCell, sim_waveA

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
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData

  real :: xcell,xedge,ycell,yedge

  real :: A0

  real :: A, B, emp, fs, x0, y0, r0

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

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,CENTER)

  A0 = sim_waveA

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

          !if (ycell .LE. 0.0) then
          !   solnData(DFUN_VAR,i,j,k) = 0.0 - ycell
          !else
          !   solnData(DFUN_VAR,i,j,k) = -1.0 * ycell
          !end if
          !solnData(DFUN_VAR,i,j,k) = ycell
          
!          solnData(DFUN_VAR,i,j,k) = ycell - A0/(cosh(sqrt(3.0*A0)*xcell/2)**2)

          ! Ellipse with disturbution
          ! phi = fs*(sqrt(x^2/A^2+y^2/B^2)-1)
          ! fs = emp + (x-x0)^2 + (y-y0)^2
          !A = 4.0d0
          !B = 2.0d0
          !emp = 0.1d0
          !x0 = 3.5d0
          !y0 = 2.0d0
          !fs = emp + (xcell-x0)**2 + (ycell-y0)**2
          !fs = 0.125*fs
          !solnData(DFUN_VAR,i,j,k) = fs*(sqrt((xcell/A)**2+(ycell/B)**2)-1.0d0)
          !solnData(CURV_VAR,i,j,k) = (sqrt((xcell/A)**2+(ycell/B)**2)-1.0d0)
          !solnData(VISC_VAR,i,j,k) = fs

          ! Circle with distrubution
          ! phi = fs*(sqrt(x^2/A^2+y^2/B^2)-1)
          ! fs = emp + (x-x0)^2 + (y-y0)^2
          !A = 3.0d0
          !B = 3.0d0
          !emp = 0.1d0
          !x0 = 1.7d0
          !y0 = 1.7d0
          !fs = emp + (xcell-x0)**2 + (ycell-y0)**2
          !fs = 1.0d0*fs
          !solnData(DFUN_VAR,i,j,k) = fs*(sqrt((xcell/A)**2+(ycell/B)**2)-1.0d0)

          ! Shizhao 
          ! Jul 7, 2015
          ! Zalesak's circle (1979, JCP)
          !r0 = 0.15d0
          !x0 = 0.0d0
          !y0 = 0.25d0
          !solnData(DFUN_VAR,i,j,k) = r0 - sqrt((xcell-x0)**2+(ycell-y0)**2)
          !if (sqrt((xcell-x0)**2+(ycell-y0)**2) < r0) then
          !  if(abs(xcell-x0) < 0.03d0 .and. (ycell-y0) < 0.1d0) then
          !    solnData(DFUN_VAR,i,j,k) = -1.0d0
          !  endif
          !endif

          ! Shizhao 
          ! Jul 7, 2015
          ! Circle
           r0 = 1.0d0 !0.5d0
           x0 = 0.0d0
           y0 = 0.0d0
           solnData(DFUN_VAR,i,j,k) = r0 - sqrt((xcell-x0)**2+(ycell-y0)**2)

        enddo
     enddo
  enddo

#if(0)
  !- wsz - Initialize the velocity in the 1st quadrant 
  do k=1,1
     do j=1,blkLimitsGC(HIGH,JAXIS)
        do i=1,blkLimitsGC(HIGH,IAXIS)+1

           ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                   real(j - NGUARD - 1)*del(JAXIS)  +  &
                   0.5*del(JAXIS)

          if (ycell .LE. 0.0) then
             facexData(VELC_FACE_VAR,i,j,k) = 1.0d0
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

  solnData(CURV_VAR,:,:,:) = 0.0
  solnData(SIGP_VAR,:,:,:) = 0.0
  solnData(VISC_VAR,:,:,:) = 0.0
  solnData(PFUN_VAR,:,:,:) = 0.0
  !solnData(DFUN_VAR,:,:,:) = 0.0

  facexData(VELC_FACE_VAR,:,:,:) = 0.0
  faceyData(VELC_FACE_VAR,:,:,:) = 0.0
  facexData(RHDS_FACE_VAR,:,:,:) = 0.0
  faceyData(RHDS_FACE_VAR,:,:,:) = 0.0

  facexData(SIGM_FACE_VAR,:,:,:) = 0.0
  faceyData(SIGM_FACE_VAR,:,:,:) = 0.0
  facexData(RH1F_FACE_VAR,:,:,:) = 0.0
  faceyData(RH1F_FACE_VAR,:,:,:) = 0.0
  facexData(RH2F_FACE_VAR,:,:,:) = 0.0
  faceyData(RH2F_FACE_VAR,:,:,:) = 0.0


!!$  ! Point to blocks center and face vars:
!!$  call Grid_getBlkPtr(blockID,solnData,CENTER)
!!$  call Grid_getBlkPtr(blockID,facexData,FACEX)
!!$  call Grid_getBlkPtr(blockID,faceyData,FACEY)
!!$
!!$
!!$  do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
!!$     do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
!!$     do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
!!$
!!$
!!$     if (ISNAN(facexData(VELC_FACE_VAR,i,j,k))) then
!!$       write(*,*) 'facexData block=',blockID
!!$       write(*,*) 'i,j,k=',i,j,k,' is a NAN.',facexData(VELC_FACE_VAR,i,j,k)
!!$     endif
!!$
!!$     if (ISNAN(faceyData(VELC_FACE_VAR,i,j,k))) then
!!$       write(*,*) 'faceyData block=',blockID
!!$       write(*,*) 'i,j,k=',i,j,k,' is a NAN.',faceyData(VELC_FACE_VAR,i,j,k)
!!$     endif
!!$
!!$
!!$     enddo
!!$     enddo
!!$  enddo

  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  call Grid_releaseBlkPtr(blockID,facexData,FACEX)
  call Grid_releaseBlkPtr(blockID,faceyData,FACEY)



  return

111    format (i4,3x,i4)
112    format (3(3x,e12.4))

end subroutine Simulation_initBlock
