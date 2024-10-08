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

  real :: A, B, emp, fs, x0, y0, r0, solnX, x1, y1, x2, y2, d1, d2, d3,&
          r_test,d_buf

  real :: x3,x4,x5,x6,y3,y4,y5,y6
  real :: d4,d5,d6,d7
  real :: fn(8)
 
  !----------------------------------------------------------------------
  
  !if (myPE .eq. MASTER_PE) write(*,*) 'InitBlockTime =',dr_simTime

  ! Get nxb, nyb and nxb:
  !call Grid_getBlkIndexSize(blockId,blIndSize,blIndSizeGC)

  !nxb = blIndSize(1)
  !nyb = blIndSize(2)
  !nzb = blIndSize(3)

  fn(1) = -1938.3
  fn(2) = -958.9
  fn(3) =  3255.2
  fn(4) = -1911.5
  fn(5) =  420.5
  fn(6) = -10.0
  fn(7) = -09.2
  fn(8) =  1.0

  !fn(1) = -2306.64
  !fn(2) =  2719.50
  !fn(3) = -866.15
  !fn(4) = -130.75
  !fn(5) =  106.30
  !fn(6) = -1.25
  !fn(7) = -6.87
  !fn(8) =  0.99

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
  solnX = 0.50007326145904204295640899471226

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

           !r0 = 0.25d0 !1.0d0 !0.5d0
           !r0 = 0.03
           !r0 = 1.0e-4
           !r0 = 0.5d0
           !r0 = 5.0e-5
           !r0 = 0.02e-3

           !r0 = 0.08
           !r0 = 0.1
           !r0 = 0.1
           r0 = 0.5
           !r_test = 0.25

           !x0 = 0.12d0
           !y0 = r0*cos((25.0/180.0)*acos(-1.0))
           y0 = 4.75d0

           !x1 = -0.12d0
           y1 = 1.5
           !y1 = r0*cos((54.0/180.0)*acos(-1.0))
           !y1 = r0*cos((35.0/180.0)*acos(-1.0))
           !y1 = r0*cos((38.0/180.0)*acos(-1.0))

           x0 =  0.0d0
           x1 =  0.0d0
           x2 =  0.0d0
           x3 =  0.9d0
           x4 = -0.3d0
           x5 = -0.6d0
           x6 = -0.9d0

           y2 = y1
           y3 = y1
           y4 = y1
           y5 = y1
           y6 = y1
           !y1 = 3.0d0


           d1 = r_test - sqrt((xcell-x0)**2+(ycell-y0)**2)
           d2 = r_test - sqrt((xcell-x1)**2+(ycell-y1)**2) 
           d3 = r0 - sqrt((xcell-x2)**2+(ycell-y2)**2)
           d4 = r0 - sqrt((xcell-x3)**2+(ycell-y3)**2)
           d5 = r0 - sqrt((xcell-x4)**2+(ycell-y4)**2)
           d6 = r0 - sqrt((xcell-x5)**2+(ycell-y5)**2)
           d7 = r0 - sqrt((xcell-x6)**2+(ycell-y6)**2)

           !d_buf = ycell-(sim_yMax-sim_yMin-3.0)
           d_buf = 2.0 - sqrt((xcell-0.0d0)**2+(ycell-9.5d0)**2)

           solnData(DFUN_VAR,i,j,k) = d3

           !if(abs(d1)<abs(d3)) then
           ! solnData(DFUN_VAR,i,j,k) = d1
           !else 
           ! solnData(DFUN_VAR,i,j,k) = d3
           !end if

           !if(abs(d1) < abs(d2) .and. abs(d1)<abs(d3)) solnData(DFUN_VAR,i,j,k) = d1
           !if(abs(d2) < abs(d3) .and. abs(d2)<abs(d1)) solnData(DFUN_VAR,i,j,k) = d2
           !if(abs(d3) < abs(d1) .and. abs(d3)<abs(d2)) solnData(DFUN_VAR,i,j,k) = d3
           
           !if(abs(d1) < abs(d2) .and. abs(d1) < abs(d3) .and. abs(d1) < abs(d4) .and. &
           !   abs(d1) < abs(d5) .and. abs(d1) < abs(d6) .and. abs(d1) < abs(d7)) solnData(DFUN_VAR,i,j,k) = d1

           !if(abs(d2) < abs(d1) .and. abs(d2) < abs(d3) .and. abs(d2) < abs(d4) .and. &
           !   abs(d2) < abs(d5) .and. abs(d2) < abs(d6) .and. abs(d2) < abs(d7)) solnData(DFUN_VAR,i,j,k) = d2

           !if(abs(d3) < abs(d2) .and. abs(d3) < abs(d1) .and. abs(d3) < abs(d4) .and. &
           !   abs(d3) < abs(d5) .and. abs(d3) < abs(d6) .and. abs(d3) < abs(d7)) solnData(DFUN_VAR,i,j,k) = d3

           !if(abs(d4) < abs(d2) .and. abs(d4) < abs(d3) .and. abs(d4) < abs(d1) .and. &
           !   abs(d4) < abs(d5) .and. abs(d4) < abs(d6) .and. abs(d4) < abs(d7)) solnData(DFUN_VAR,i,j,k) = d4

           !if(abs(d5) < abs(d2) .and. abs(d5) < abs(d3) .and. abs(d5) < abs(d4) .and. &
           !   abs(d5) < abs(d1) .and. abs(d5) < abs(d6) .and. abs(d5) < abs(d7)) solnData(DFUN_VAR,i,j,k) = d5

           !if(abs(d6) < abs(d2) .and. abs(d6) < abs(d3) .and. abs(d6) < abs(d4) .and. &
           !   abs(d6) < abs(d5) .and. abs(d6) < abs(d1) .and. abs(d6) < abs(d7)) solnData(DFUN_VAR,i,j,k) = d6

           !if(abs(d7) < abs(d2) .and. abs(d7) < abs(d3) .and. abs(d7) < abs(d4) .and. &
           !   abs(d7) < abs(d5) .and. abs(d7) < abs(d6) .and. abs(d7) < abs(d1)) solnData(DFUN_VAR,i,j,k) = d7

           !solnData(TEMP_VAR,i,j,k) = 0.1185 + (-0.1185/erf(solnX))*(erf(ycell)/(2*sqrt(0.25)))
           !solnData(DFUN_VAR,i,j,k) = sqrt((xcell-x0)**2+(ycell-y0)**2) - r0
           !solnData(DFUN_VAR,i,j,k) = (0.08/128. )*(4 + cos(2*acos(-1.0)*(xcell-(0.08/2.))/0.08)) - ycell
           !solnData(DFUN_VAR,i,j,k) = (1.0/128. )*(4 + cos(2*acos(-1.0)*(xcell-(1.0/2.0))/1.0)) - ycell
           !solnData(DFUN_VAR,i,j,k) = (0.0273/128. )*(4 + cos(2*acos(-1.0)*(xcell-(0.0273/2.))/0.0273)) - ycell
           !solnData(DFUN_VAR,i,j,k) = (0.0023/128.0)*(4 + cos(2*acos(-1.0)*(xcell-(0.0023/2.))/0.0023)) - ycell
           !solnData(DFUN_VAR,i,j,k) = (0.08/128. )*(4 + cos(2*acos(-1.0)*(xcell)/0.08)) - ycell
           !solnData(DFUN_VAR,i,j,k) = 0.5 - ycell

           solnData(TEMP_VAR,i,j,k) = 0.0

           !if(ycell .le. 0.3520 .and. solnData(DFUN_VAR,i,j,k) .lt. 0.0) solnData(TEMP_VAR,i,j,k) = (0.3520-ycell)/0.3520
         
           !if(ycell .le. 0.2 .and. solnData(DFUN_VAR,i,j,k) .lt. 0.0) solnData(TEMP_VAR,i,j,k) = (0.2-ycell)/0.2

           !if(ycell .le. 0.5 .and. solnData(DFUN_VAR,i,j,k) .lt. 0.0) solnData(TEMP_VAR,i,j,k) = (0.5-ycell)/0.5
           !if(solnData(DFUN_VAR,i,j,k) .ge. 0.0) solnData(TEMP_VAR,i,j,k) = 0.4

           !if(ycell .le. 0.3781 .and. solnData(DFUN_VAR,i,j,k) .lt. 0.0) solnData(TEMP_VAR,i,j,k) = (0.3781-ycell)/0.3781

           !if(ycell .le. 0.3520 .and. solnData(DFUN_VAR,i,j,k) .lt. 0.0) then

           !solnData(TEMP_VAR,i,j,k) = fn(1)*(ycell**7) + fn(2)*(ycell**6) + fn(3)*(ycell**5) + &
           !                           fn(4)*(ycell**4) + fn(5)*(ycell**3) + fn(6)*(ycell**2) + &
           !                           fn(7)*(ycell**1) + fn(8)

           !if (solnData(TEMP_VAR,i,j,k) .lt. 0.0) solnData(TEMP_VAR,i,j,k) = 0.0


           !end if

           !if(ycell .le. 0.3792 .and. solnData(DFUN_VAR,i,j,k) .lt. 0.0) then

           !solnData(TEMP_VAR,i,j,k) = fn(1)*(ycell**7) + fn(2)*(ycell**6) + fn(3)*(ycell**5) + &
           !                           fn(4)*(ycell**4) + fn(5)*(ycell**3) + fn(6)*(ycell**2) + &
           !                           fn(7)*(ycell**1) + fn(8)

           !if (solnData(TEMP_VAR,i,j,k) .lt. 0.0) solnData(TEMP_VAR,i,j,k) = 0.0


           !end if

           !if(solnData(DFUN_VAR,i,j,k) .ge. 0.0) solnData(TEMP_VAR,i,j,k) = 0.5

           !if(ycell .le. 0.3792 .and. solnData(DFUN_VAR,i,j,k) .lt. 0.0) solnData(TEMP_VAR,i,j,k) = (0.3792-ycell)/0.3792

           !if(ycell .le. 0.2 .and. solnData(DFUN_VAR,i,j,k) .lt. 0.0) solnData(TEMP_VAR,i,j,k) = (0.2-ycell)/0.2

           !if(ycell .le. 0.0044) solnData(TEMP_VAR,i,j,k) = (0.0044-ycell)/0.0044
 
           !if(ycell .le. 8.7308) solnData(TEMP_VAR,i,j,k) = (8.7308 - ycell)/8.7308
     
           !if(ycell .le. 9.9714 .and. solnData(DFUN_VAR,i,j,k) .lt. 0.0) solnData(TEMP_VAR,i,j,k) = (9.9714 - ycell)/9.9714

           !if(ycell .le. 9.9035 .and. solnData(DFUN_VAR,i,j,k) .lt. 0.0) solnData(TEMP_VAR,i,j,k) = (9.9035 - ycell)/9.9035
           !if(solnData(DFUN_VAR,i,j,k) .ge. 0.0) solnData(TEMP_VAR,i,j,k) = 0.3578

           !if(ycell .le. 9.3178 .and. solnData(DFUN_VAR,i,j,k) .lt. 0.0) solnData(TEMP_VAR,i,j,k) = (9.3178 - ycell)/9.3178

           !if(solnData(DFUN_VAR,i,j,k) .lt. 0.0) then

           !  if(ycell .le. 0.0800) solnData(TEMP_VAR,i,j,k) = (2.0 - ycell)/2.0

           !end if

           !if(solnData(DFUN_VAR,i,j,k) .ge. 0.) then
             !solnData(TEMP_VAR,i,j,k) = 0.1*(0.08-ycell)/0.08
           !  solnData(TEMP_VAR,i,j,k) = (solnData(DFUN_VAR,i,j,k)*1.0)/(solnData(DFUN_VAR,i,j,k)+ycell)
           !else
           !  solnData(TEMP_VAR,i,j,k) = 0.0
           !end if

           !solnData(TEMP_VAR,i,j,k) = 0.0          
 
           !solnData(DFUN_VAR,i,j,k) = 0.5 - ycell

           !if(solnData(DFUN_VAR,i,j,k) .ge. 0.) then
           !  solnData(TEMP_VAR,i,j,k) = 0.1185 + (-0.1185/erf(solnX))*(erf(ycell)/(2*sqrt(0.25)))
           !else
           !  solnData(TEMP_VAR,i,j,k) = 0.0
           !end if

           !if (solnData(DFUN_VAR,i,j,k) .ge. 0) then

              !solnData(TEMP_VAR,i,j,k) = 0.0
           !   faceyData(VELC_FACE_VAR,i,j,k) = 0.0

           !else

           !   solnData(TEMP_VAR,i,j,k) = 0.0
           !   faceyData(VELC_FACE_VAR,i,j,k) = 0.0

           !end if

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

  solnData(CURV_VAR,:,:,:) = 0.0
  solnData(SIGP_VAR,:,:,:) = 0.0
  solnData(VISC_VAR,:,:,:) = 0.0
  solnData(PFUN_VAR,:,:,:) = 0.0
  !solnData(DFUN_VAR,:,:,:) = 0.0
  solnData(NRMX_VAR,:,:,:) = 0.0
  solnData(NRMY_VAR,:,:,:) = 0.0
  solnData(NRMZ_VAR,:,:,:) = 0.0
  solnData(SMHV_VAR,:,:,:) = 0.0
  solnData(SMRH_VAR,:,:,:) = 0.0
  solnData(MDOT_VAR,:,:,:) = 0.0
  solnData(RHST_VAR,:,:,:) = 0.0
  solnData(TOLD_VAR,:,:,:) = 0.0

  solnData(PTES_VAR,:,:,:) = 0.0
  solnData(RTES_VAR,:,:,:) = 0.0
  solnData(ALPH_VAR,:,:,:) = 0.0

  facexData(VELC_FACE_VAR,:,:,:) = 0.0
  faceyData(VELC_FACE_VAR,:,:,:) = 0.0
  facexData(VELI_FACE_VAR,:,:,:) = 0.0
  faceyData(VELI_FACE_VAR,:,:,:) = 0.0
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
