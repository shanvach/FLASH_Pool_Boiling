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

  x0 = 0.0
  y0 = 0.0
  r0 = 0.5

  solnData(DELE_VAR,:,:,:)  = 1e6          ! particles/m3
  solnData(PRHV_VAR,:,:,:)  = 101325.0     ! Pascal, heavy particle pressure  
  solnData(TPHV_VAR,:,:,:)  = 300.0        ! Kelvin, heavy particle temperature
  solnData(TPEL_VAR,:,:,:)  = 0.8*11604.52 ! Kelvin, electron temperature
  solnData(DNAT_VAR,:,:,:)  = 0.0          ! total neutral species
  solnData(DNIT_VAR,:,:,:)  = 0.0          ! total ion species
  solnData(GNE_VAR,:,:,:)   = 0.0          ! electron generation rate
  solnData(GNEBZ_VAR,:,:,:) = 0.0          ! electron generation rate (Boltzmann)
  solnData(GNERT_VAR,:,:,:) = 0.0          ! ratio
 
  !initialize values of heavy species and generation rates
  do i=0,9
    solnData(DHV0_VAR+i,:,:,:) = 1e6
    solnData(GNH0_VAR+i,:,:,:) = 0.0
  end do
 
  !initialize values for reaction rates of heavy species
  do i=0,16
    solnData(RSP0_VAR+i,:,:,:) = 0.0
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

           
           solnData(DFUN_VAR,i,j,k) = r0 - sqrt((xcell-x0)**2+(ycell-y0)**2)

           if(solnData(DFUN_VAR,i,j,k) .ge. 0.0) then

                solnData(DELE_VAR,i,j,k) = 1e18 ! Electrons
                solnData(DHV0_VAR,i,j,k) = 0.9*1e26 ! He
                solnData(DHV6_VAR,i,j,k) = 1e18 ! He+
                solnData(DHV1_VAR,i,j,k) = 0.1*0.8*1e26 ! N2 
                solnData(DHV2_VAR,i,j,k) = 0.1*0.8*1e26 ! O2

           end if

          enddo
     enddo
  enddo
  !
  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  call Grid_releaseBlkPtr(blockID,facexData,FACEX)
  call Grid_releaseBlkPtr(blockID,faceyData,FACEY)

  return

111    format (i4,3x,i4)
112    format (3(3x,e12.4))

end subroutine Simulation_initBlock
