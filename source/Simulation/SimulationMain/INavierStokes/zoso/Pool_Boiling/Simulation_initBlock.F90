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
                              sim_nuc_radii

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

  real :: A, B, emp, fs, x0, y0, r0, solnX, x1, y1, x2, y2, d1, d2, d3, r_test
  real :: x3,x4,x5,x6,y3,y4,y5,y6
  real :: d4,d5,d6,d7
  real :: nuc_dfun
  integer :: nuc_index
  real :: fn(8)
 
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
  solnX = 0.50007326145904204295640899471226

  sim_nuc_radii(1) = 0.12
  sim_nuc_radii(2) = 0.17
  sim_nuc_radii(3) = 0.15
  sim_nuc_radii(4) = 0.10
  sim_nuc_radii(5) = 0.30
  sim_nuc_radii(6) = 0.12
  sim_nuc_radii(7) = 0.15
  sim_nuc_radii(8) = 0.17
  sim_nuc_radii(9) = 0.10


  sim_nuc_site_x(1) = -3.00
  sim_nuc_site_x(2) = -2.25
  sim_nuc_site_x(3) = -1.50
  sim_nuc_site_x(4) = -0.75
  sim_nuc_site_x(5) =  0.00
  sim_nuc_site_x(6) =  0.75
  sim_nuc_site_x(7) =  1.50
  sim_nuc_site_x(8) =  2.25
  sim_nuc_site_x(9) =  3.00

  do nuc_index=1,9

        sim_nuc_site_y(nuc_index) = sim_nuc_radii(nuc_index)*cos((30.0/180.0)*acos(-1.0))

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

           do nuc_index=1,9

           nuc_dfun  = sim_nuc_radii(nuc_index) - sqrt((xcell-sim_nuc_site_x(nuc_index))**2+(ycell-sim_nuc_site_y(nuc_index))**2);

           if (nuc_index == 1) then

                solnData(DFUN_VAR,i,j,k) = nuc_dfun

           else

                solnData(DFUN_VAR,i,j,k) = max(solnData(DFUN_VAR,i,j,k),nuc_dfun)

           end if

           end do

           solnData(TEMP_VAR,i,j,k) = 0.0

           !if(ycell .le. 2.0 .and. solnData(DFUN_VAR,i,j,k) .lt. 0.0) solnData(TEMP_VAR,i,j,k) = (2.0 - ycell)/2.0
           if(ycell .le. 2.0 .and. xcell .ge. -4.0 .and. xcell .le. 4.0) solnData(TEMP_VAR,i,j,k) = (2.0 - ycell)/2.0

          enddo
     enddo
  enddo

  sim_nuc_site_y(:) = 0.1*cos((30.0/180.0)*acos(-1.0))

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
