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

  real :: xcell,xedge,ycell,yedge,zcell

  real :: A0

  real :: A,B,emp,fs,x0,y0,r0,solnX,z0,x1,y1,z1,x2,y2,z2,d1,d2,d3

  real :: x3,y3,z3,d4,r1

  real :: fn(8)
  real :: x4,x5,x6,x7,x8
  real :: y4,y5,y6,y7,y8
  real :: d5,d6,d7,d8,d9
  real :: z4,z5,z6,z7,z8 

  real, dimension(17) :: Nuc_radii,Nuc_sites_x,Nuc_sites_z,Nuc_sites_y
  real :: Nuc_dfun
  integer :: Nuc_Index

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

  Nuc_radii(1) =  0.0213
  Nuc_radii(2) =  0.1012
  Nuc_radii(3) =  0.0438
  Nuc_radii(4) =  0.0258
  Nuc_radii(5) =  0.0225
  Nuc_radii(6) =  0.0516
  Nuc_radii(7) =  0.0236
  Nuc_radii(8) =  0.0505
  Nuc_radii(9) =  0.1245
  Nuc_radii(10) =  0.0236
  Nuc_radii(11) =  0.0247
  Nuc_radii(12) =  0.0203
  Nuc_radii(13) =  0.0248
  Nuc_radii(14) =  0.0303
  Nuc_radii(15) =  0.0236
  Nuc_radii(16) =  0.1301
  Nuc_radii(17) =  0.1021

  
  Nuc_sites_x(1) = -0.7844
  Nuc_sites_x(2) = -0.6301
  Nuc_sites_x(3) = -0.5038
  Nuc_sites_x(4) = -0.2685
  Nuc_sites_x(5) = -0.2426
  Nuc_sites_x(6) = -0.1658
  Nuc_sites_x(7) = -0.1304
  Nuc_sites_x(8) =  0.0631
  Nuc_sites_x(9) =  0.1092
  Nuc_sites_x(10) =  0.1645
  Nuc_sites_x(11) =  0.1715
  Nuc_sites_x(12) =  0.1790
  Nuc_sites_x(13) =  0.2413
  Nuc_sites_x(14) =  0.3836
  Nuc_sites_x(15) =  0.4113
  Nuc_sites_x(16) =  0.5302
  Nuc_sites_x(17) =  0.5342


  Nuc_sites_z(1) =  0.3922
  Nuc_sites_z(2) =  -0.3964
  Nuc_sites_z(3) =  -0.2029
  Nuc_sites_z(4) =   -0.2934
  Nuc_sites_z(5) =  -0.5689
  Nuc_sites_z(6) =  -0.4175
  Nuc_sites_z(7) =  -0.5689
  Nuc_sites_z(8) =  -0.4027
  Nuc_sites_z(9) =   0.6362
  Nuc_sites_z(10) =  -0.2787
  Nuc_sites_z(11) =  -0.1756
  Nuc_sites_z(12) =   0.0410
  Nuc_sites_z(13) =  -0.1230
  Nuc_sites_z(14) =  -0.4763
  Nuc_sites_z(15) =  -0.2639
  Nuc_sites_z(16) =   0.4175
  Nuc_sites_z(17) =  -0.4111

  do Nuc_Index=1,17

       Nuc_sites_y(Nuc_Index) = Nuc_radii(Nuc_Index)*cos((35.0/180.0)*acos(-1.0))

  end do

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

           zcell  = coord(KAXIS) - bsize(KAXIS)/2.0 +  &
                   real(k - NGUARD - 1)*del(KAXIS)  +  &
                   0.5*del(KAXIS)


           !_________________TEST PROBLEM 1_____________________!

           !! Main Bubble
           !!r0 =  0.05
           !r0 =  0.1
           r0 =  3.3599
           !r0 = 0.08
           x0 =  0.0
           z0 =  0.0
           !y0 =  r0*cos((30.0/180.0)*acos(-1.0))
           !!y0 =  r0*cos((54.0/180.0)*acos(-1.0))
           y0 =  r0*cos((35.0/180.0)*acos(-1.0))

           !! Auxiallary Bubbles
           !r1 =  0.1
          
           !!y1 = r1*cos((30.0/180.0)*acos(-1.0))
           !y1 = r1*cos((35.0/180.0)*acos(-1.0))
           !y2 = y1
           !y3 = y1
           !y4 = y1
           !y5 = y1
           !y6 = y1
           !y7 = y1
           !y8 = y1

           !!x1 =   0.0
           !!x2 =   0.3/sqrt(2.)
           !!x3 =   0.3
           !!x4 =   0.3/sqrt(2.)
           !!x5 =   0.0
           !!x6 =  -0.3/sqrt(2.)
           !!x7 =  -0.3
           !!x8 =  -0.3/sqrt(2.)

           !!z1 = -sqrt(0.09 - x1**2)
           !!z2 = -sqrt(0.09 - x2**2)
           !!z3 =  sqrt(0.09 - x3**2)
           !!z4 =  sqrt(0.09 - x4**2)
           !!z5 =  sqrt(0.09 - x5**2)
           !!z6 =  sqrt(0.09 - x6**2)
           !!z7 =  sqrt(0.09 - x7**2)
           !!z8 = -sqrt(0.09 - x8**2)

           !x0 =  3.5183
           !x1 = -0.7925
           !x2 =  1.7887
           !x3 = -2.7197
           !x4 = -4.2251
           !x5 =  1.1967
           !x6 = -4.0002
           !x7 = -0.1807
           !x8 =  3.8425

           !z0 =  0.9888
           !z1 = -0.9090
           !z2 =  3.2350
           !z3 =  2.7494
           !z4 =  0.6905
           !z5 = -2.8537
           !z6 = -2.3406
           !z7 =  3.4786
           !z8 = -2.2419

           !! Distance functions
           d1 = r0 - sqrt((xcell-x0)**2+(ycell-y0)**2+(zcell-z0)**2)
           !d2 = r1 - sqrt((xcell-x1)**2+(ycell-y1)**2+(zcell-z1)**2)
           !d3 = r1 - sqrt((xcell-x2)**2+(ycell-y2)**2+(zcell-z2)**2)
           !d4 = r1 - sqrt((xcell-x3)**2+(ycell-y3)**2+(zcell-z3)**2)
           !d5 = r1 - sqrt((xcell-x4)**2+(ycell-y4)**2+(zcell-z4)**2)
           !d6 = r1 - sqrt((xcell-x5)**2+(ycell-y5)**2+(zcell-z5)**2)
           !d7 = r1 - sqrt((xcell-x6)**2+(ycell-y6)**2+(zcell-z6)**2)
           !d8 = r1 - sqrt((xcell-x7)**2+(ycell-y7)**2+(zcell-z7)**2)
           !d9 = r1 - sqrt((xcell-x8)**2+(ycell-y8)**2+(zcell-z8)**2)

           !! Single bubble setup
           solnData(DFUN_VAR,i,j,k) = d1

           !! Multiple bubble setup         
           
           !if(abs(d1) < abs(d2) .and. abs(d1) < abs(d3) .and. abs(d1) < abs(d4) .and. &
           !   abs(d1) < abs(d5) .and. abs(d1) < abs(d6) .and. abs(d1) < abs(d7) .and. &
           !   abs(d1) < abs(d8) .and. abs(d1) < abs(d9)) solnData(DFUN_VAR,i,j,k) = d1

           !if(abs(d2) < abs(d1) .and. abs(d2) < abs(d3) .and. abs(d2) < abs(d4) .and. &
           !   abs(d2) < abs(d5) .and. abs(d2) < abs(d6) .and. abs(d2) < abs(d7) .and. &
           !   abs(d2) < abs(d8) .and. abs(d2) < abs(d9)) solnData(DFUN_VAR,i,j,k) = d2

           !if(abs(d3) < abs(d2) .and. abs(d3) < abs(d1) .and. abs(d3) < abs(d4) .and. &
           !   abs(d3) < abs(d5) .and. abs(d3) < abs(d6) .and. abs(d3) < abs(d7) .and. &
           !   abs(d3) < abs(d8) .and. abs(d3) < abs(d9)) solnData(DFUN_VAR,i,j,k) = d3

           !if(abs(d4) < abs(d2) .and. abs(d4) < abs(d3) .and. abs(d4) < abs(d1) .and. &
           !   abs(d4) < abs(d5) .and. abs(d4) < abs(d6) .and. abs(d4) < abs(d7) .and. &
           !   abs(d4) < abs(d8) .and. abs(d4) < abs(d9)) solnData(DFUN_VAR,i,j,k) = d4

           !if(abs(d5) < abs(d2) .and. abs(d5) < abs(d3) .and. abs(d5) < abs(d4) .and. &
           !   abs(d5) < abs(d1) .and. abs(d5) < abs(d6) .and. abs(d5) < abs(d7) .and. &
           !   abs(d5) < abs(d8) .and. abs(d5) < abs(d9)) solnData(DFUN_VAR,i,j,k) = d5

           !if(abs(d6) < abs(d2) .and. abs(d6) < abs(d3) .and. abs(d6) < abs(d4) .and. &
           !   abs(d6) < abs(d5) .and. abs(d6) < abs(d1) .and. abs(d6) < abs(d7) .and. &
           !   abs(d6) < abs(d8) .and. abs(d6) < abs(d9)) solnData(DFUN_VAR,i,j,k) = d6

           !if(abs(d7) < abs(d2) .and. abs(d7) < abs(d3) .and. abs(d7) < abs(d4) .and. &
           !   abs(d7) < abs(d5) .and. abs(d7) < abs(d6) .and. abs(d7) < abs(d1) .and. &
           !   abs(d7) < abs(d8) .and. abs(d7) < abs(d9)) solnData(DFUN_VAR,i,j,k) = d7

           !if(abs(d8) < abs(d2) .and. abs(d8) < abs(d3) .and. abs(d8) < abs(d4) .and. &
           !   abs(d8) < abs(d5) .and. abs(d8) < abs(d6) .and. abs(d8) < abs(d7) .and. &
           !   abs(d8) < abs(d1) .and. abs(d8) < abs(d9)) solnData(DFUN_VAR,i,j,k) = d8

           !if(abs(d9) < abs(d2) .and. abs(d9) < abs(d3) .and. abs(d9) < abs(d4) .and. &
           !   abs(d9) < abs(d5) .and. abs(d9) < abs(d6) .and. abs(d9) < abs(d7) .and. &
           !   abs(d9) < abs(d8) .and. abs(d9) < abs(d1)) solnData(DFUN_VAR,i,j,k) = d9

           solnData(TEMP_VAR,i,j,k) = 0.0

           if(ycell .le. 9.7721 .and. solnData(DFUN_VAR,i,j,k) .lt. 0.0) solnData(TEMP_VAR,i,j,k) = (9.7721 - ycell)/9.7721
           !if(ycell .le. 0.3520 .and. solnData(DFUN_VAR,i,j,k) .lt. 0.0) solnData(TEMP_VAR,i,j,k) = (0.3520 - ycell)/0.3520

           !if(solnData(DFUN_VAR,i,j,k) .ge. 0.0) solnData(TEMP_VAR,i,j,k) = 0.1

           !if(ycell .le. 0.3520 .and. solnData(DFUN_VAR,i,j,k) .lt. 0.0) then
           !!!if(ycell .le. 0.3792 .and. solnData(DFUN_VAR,i,j,k) .lt. 0.0) then

           !solnData(TEMP_VAR,i,j,k) = fn(1)*(ycell**7) + fn(2)*(ycell**6) + fn(3)*(ycell**5) + &
           !                           fn(4)*(ycell**4) + fn(5)*(ycell**3) + fn(6)*(ycell**2) + &
           !                           fn(7)*(ycell**1) + fn(8)

           !if (solnData(TEMP_VAR,i,j,k) .lt. 0.0) solnData(TEMP_VAR,i,j,k) = 0.0


           !end if

           !if(solnData(DFUN_VAR,i,j,k) .ge. 0.0) solnData(TEMP_VAR,i,j,k) = 0.5

           !!_______________PRODUCTION RUN PROBLEM 1_____________________!

           !solnData(DFUN_VAR,i,j,k) = 1E10

           !do Nuc_Index=1,17

           ! Nuc_dfun  = Nuc_radii(Nuc_Index) -  sqrt((xcell-Nuc_sites_x(Nuc_Index))**2+(ycell-Nuc_sites_y(Nuc_Index))**2+(zcell-Nuc_sites_z(Nuc_Index))**2);
            
           ! if(abs(solnData(DFUN_VAR,i,j,k)) > abs(Nuc_dfun)) solnData(DFUN_VAR,i,j,k) = Nuc_dfun
                
           !end do

           !solnData(TEMP_VAR,i,j,k) = -2.0769
           !if(ycell .le. 0.4 .and. solnData(DFUN_VAR,i,j,k) .lt. 0.0) solnData(TEMP_VAR,i,j,k) = (((0.4 - ycell)*1.0)+(ycell*(-2.0769)))/0.4
           !if(solnData(DFUN_VAR,i,j,k) .ge. 0.0) solnData(TEMP_VAR,i,j,k) = 0.0
        

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
