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

  use Simulation_Data, Only : sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_init

  use Grid_interface, only  : Grid_getCellMetrics, Grid_getBlkIndexLimits, &
                              Grid_getCellCoords, &
                              Grid_getBlkCenterCoords, Grid_getBlkBoundBox, &
                              Grid_getBlkPtr, Grid_releaseBlkPtr

  use Driver_Data, Only     : dr_simTime

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: blockID
 
  real, pointer, dimension(:,:,:,:) :: solnData, facexData, faceyData, facezData
  
  integer :: i, j, k
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM)   :: blkIndSize, blkIndSizeGC
  real, dimension(MDIM)      :: coord, bsize
  real, dimension(2,MDIM)    :: boundBox
  real                       :: pi
  real, dimension(GRID_IHI_GC) :: dx
  real, dimension(GRID_JHI_GC) :: dy
  real, dimension(GRID_KHI_GC) :: dz
  real, dimension(GRID_IHI_GC) :: xcell, xedge
  real, dimension(GRID_JHI_GC) :: ycell, yedge
  real, dimension(GRID_KHI_GC) :: zcell, zedge

  ! write out simulation time at initialization 
  ! write(*,*) 'BlockId =',blockId,' sim_init =',sim_init

  ! Point to Blocks centered variables:
  call Grid_getBlkPtr(blockID, solnData, CENTER)

  ! Point to Blocks face variables: 
  call Grid_getBlkPtr(blockID, facexData, FACEX)
  call Grid_getBlkPtr(blockID, faceyData, FACEY)

#if NDIM == 3
  call Grid_getBlkPtr(blockID, facezData, FACEZ)
#endif


  !Get Coord and Bsize for the block
  call Grid_getBlkBoundBox(blockId, boundBox)
  bsize(:) = boundBox(2,:) - boundBox(1,:)
  call Grid_getBlkCenterCoords(blockId, coord)

  call Grid_getCellMetrics(IAXIS, blockId, CENTER, .true., dx, GRID_IHI_GC)  
  call Grid_getCellMetrics(JAXIS, blockId, CENTER, .true., dy, GRID_JHI_GC)  
  call Grid_getCellMetrics(KAXIS, blockId, CENTER, .true., dz, GRID_KHI_GC)  

  call Grid_getCellCoords(IAXIS, blockId, CENTER, .true., xcell, GRID_IHI_GC)  
  call Grid_getCellCoords(JAXIS, blockId, CENTER, .true., ycell, GRID_JHI_GC)  
  call Grid_getCellCoords(KAXIS, blockId, CENTER, .true., zcell, GRID_KHI_GC)  
  
  call Grid_getCellCoords(IAXIS, blockId, LEFT_EDGE, .true., xedge, GRID_IHI_GC)  
  call Grid_getCellCoords(JAXIS, blockId, LEFT_EDGE, .true., yedge, GRID_JHI_GC)  
  call Grid_getCellCoords(KAXIS, blockId, LEFT_EDGE, .true., zedge, GRID_KHI_GC)  
  

  ! Initialize Simulation Blocks
  select case (sim_init)



    
    ! Rayleigh Benard flow (Cold on Hot)
    case (0)

      solnData(PRES_VAR,:,:,:) = 0.0
      solnData(DELP_VAR,:,:,:) = 0.0
      solnData(DUST_VAR,:,:,:) = 0.0
      solnData(TVIS_VAR,:,:,:) = 0.0
      solnData(TEMP_VAR,:,:,:) = 1.0

      facexData(VELC_FACE_VAR,:,:,:) = 0.0
      faceyData(VELC_FACE_VAR,:,:,:) = 0.0
      facexData(RHDS_FACE_VAR,:,:,:) = 0.0
      faceyData(RHDS_FACE_VAR,:,:,:) = 0.0

#if NDIM == 3
      facezData(VELC_FACE_VAR,:,:,:) = 0.0
      facezData(RHDS_FACE_VAR,:,:,:) = 0.0

      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, CENTER)
      do k=1, blkLimitsGC(HIGH,KAXIS)
        do j=1, blkLimitsGC(HIGH,JAXIS)
          do i=1, blkLimitsGC(HIGH,IAXIS)
            if (zcell(k) >= 0.5) then
              solnData(TEMP_VAR,i,j,k) = 0.0
            endif
          enddo
        enddo
      enddo
#endif
#if NDIM == 2
      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, CENTER)
      do k=1, blkLimitsGC(HIGH,KAXIS)
        do j=1, blkLimitsGC(HIGH,JAXIS)
          do i=1, blkLimitsGC(HIGH,IAXIS)
            if (ycell(j) >= 0.5) then
              solnData(TEMP_VAR,i,j,k) = 0.0
            endif
          enddo
        enddo
      enddo
#endif




    ! Rayleigh Benard flow (blank)
    case (1)

      solnData(PRES_VAR,:,:,:) = 0.0
      solnData(DELP_VAR,:,:,:) = 0.0
      solnData(DUST_VAR,:,:,:) = 0.0
      solnData(TVIS_VAR,:,:,:) = 0.0
      solnData(TEMP_VAR,:,:,:) = 0.5

      facexData(VELC_FACE_VAR,:,:,:) = 0.0
      faceyData(VELC_FACE_VAR,:,:,:) = 0.0
      facexData(RHDS_FACE_VAR,:,:,:) = 0.0
      faceyData(RHDS_FACE_VAR,:,:,:) = 0.0

#if NDIM == 3
      facezData(VELC_FACE_VAR,:,:,:) = 0.0
      facezData(RHDS_FACE_VAR,:,:,:) = 0.0
#endif




    ! Rayleigh Benard flow (cond)
    case (2)

      solnData(PRES_VAR,:,:,:) = 0.0
      solnData(DELP_VAR,:,:,:) = 0.0
      solnData(DUST_VAR,:,:,:) = 0.0
      solnData(TVIS_VAR,:,:,:) = 0.0

      facexData(VELC_FACE_VAR,:,:,:) = 0.0
      faceyData(VELC_FACE_VAR,:,:,:) = 0.0
      facexData(RHDS_FACE_VAR,:,:,:) = 0.0
      faceyData(RHDS_FACE_VAR,:,:,:) = 0.0

#if NDIM == 3
      facezData(VELC_FACE_VAR,:,:,:) = 0.0
      facezData(RHDS_FACE_VAR,:,:,:) = 0.0
      
      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, CENTER)
      do k=1, blkLimitsGC(HIGH,KAXIS)
        do j=1, blkLimitsGC(HIGH,JAXIS)
          do i=1, blkLimitsGC(HIGH,IAXIS)
            xcell = coord(IAXIS) - bsize(IAXIS)/2.0 + real(i-NGUARD-1)*dx(i) + 0.5*dx(i)
            zcell = coord(KAXIS) - bsize(KAXIS)/2.0 + real(k-NGUARD-1)*dz(k) + 0.5*dz(k)
!            solnData(TEMP_VAR,i,j,k) = (1.0 - zcell) 
          enddo
        enddo
      enddo
#endif
#if NDIM == 2
      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, CENTER)
      do k=1, blkLimitsGC(HIGH,KAXIS)
        do j=1, blkLimitsGC(HIGH,JAXIS)
          do i=1, blkLimitsGC(HIGH,IAXIS)
            xcell = coord(IAXIS) - bsize(IAXIS)/2.0 + real(i-NGUARD-1)*dx(i) + 0.5*dx(i)
            ycell = coord(JAXIS) - bsize(JAXIS)/2.0 + real(j-NGUARD-1)*dy(j) + 0.5*dy(j)
!            solnData(TEMP_VAR,i,j,k) = (1.0 - ycell) 
          enddo
        enddo
      enddo
#endif




    ! channel flow (left to right w/ T=0)
    case (10)

      solnData(PRES_VAR,:,:,:) = 0.0
      solnData(DELP_VAR,:,:,:) = 0.0
      solnData(DUST_VAR,:,:,:) = 0.0
      solnData(TVIS_VAR,:,:,:) = 0.0
      solnData(TEMP_VAR,:,:,:) = 0.0

      facexData(VELC_FACE_VAR,:,:,:) = 1.0
      faceyData(VELC_FACE_VAR,:,:,:) = 0.0
      facexData(RHDS_FACE_VAR,:,:,:) = 0.0
      faceyData(RHDS_FACE_VAR,:,:,:) = 0.0

#if NDIM == 3
      facezData(VELC_FACE_VAR,:,:,:) = 0.0
      facezData(RHDS_FACE_VAR,:,:,:) = 0.0
#endif      



    ! channel flow (left to right w/ init T)
    case (11)

      solnData(PRES_VAR,:,:,:) = 0.0
      solnData(DELP_VAR,:,:,:) = 0.0
      solnData(DUST_VAR,:,:,:) = 0.0
      solnData(TVIS_VAR,:,:,:) = 0.0
      solnData(TEMP_VAR,:,:,:) = 0.0

      facexData(VELC_FACE_VAR,:,:,:) = 1.0
      faceyData(VELC_FACE_VAR,:,:,:) = 0.0
      facexData(RHDS_FACE_VAR,:,:,:) = 0.0
      faceyData(RHDS_FACE_VAR,:,:,:) = 0.0

#if NDIM == 3
      facezData(VELC_FACE_VAR,:,:,:) = 0.0
      facezData(RHDS_FACE_VAR,:,:,:) = 0.0
#endif      
      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, CENTER)
      do k=1, blkLimitsGC(HIGH,KAXIS)
        do j=1, blkLimitsGC(HIGH,JAXIS)
          do i=1, blkLimitsGC(HIGH,IAXIS)
            solnData(TEMP_VAR,i,j,k) = (1.0 - xcell(i))
          enddo
        enddo
      enddo




    ! channel flow (bottom to top w/ T=8)
    case (12)

      solnData(PRES_VAR,:,:,:) = 0.0
      solnData(DELP_VAR,:,:,:) = 0.0
      solnData(DUST_VAR,:,:,:) = 0.0
      solnData(TVIS_VAR,:,:,:) = 0.0
      solnData(TEMP_VAR,:,:,:) = 8.0

      facexData(VELC_FACE_VAR,:,:,:) = 0.0
      faceyData(VELC_FACE_VAR,:,:,:) = 1.0
      facexData(RHDS_FACE_VAR,:,:,:) = 0.0
      faceyData(RHDS_FACE_VAR,:,:,:) = 0.0

#if NDIM == 3
      facezData(VELC_FACE_VAR,:,:,:) = 0.0
      facezData(RHDS_FACE_VAR,:,:,:) = 0.0
#endif      



    ! channel flow (front to back w/ T=0)
    case (13)

      solnData(PRES_VAR,:,:,:) = 0.0
      solnData(DELP_VAR,:,:,:) = 0.0
      solnData(DUST_VAR,:,:,:) = 0.0
      solnData(TVIS_VAR,:,:,:) = 0.0
      solnData(TEMP_VAR,:,:,:) = 0.0

      facexData(VELC_FACE_VAR,:,:,:) = 0.0
      faceyData(VELC_FACE_VAR,:,:,:) = 0.0
      facexData(RHDS_FACE_VAR,:,:,:) = 0.0
      faceyData(RHDS_FACE_VAR,:,:,:) = 0.0

#if NDIM == 3
      facezData(VELC_FACE_VAR,:,:,:) = 1.0
      facezData(RHDS_FACE_VAR,:,:,:) = 0.0
#endif      



    ! taylor-green flow
    case (20)
      
      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, FACEX)
      do k=1, blkLimitsGC(HIGH,KAXIS)
        do j=1, blkLimitsGC(HIGH,JAXIS)
          do i=1, blkLimitsGC(HIGH,IAXIS)
!            xedge = coord(IAXIS) - bsize(IAXIS)/2.0 + real(i-NGUARD-1)*dx(i)
!            ycell = coord(JAXIS) - bsize(JAXIS)/2.0 + real(j-NGUARD-1)*dy(j) + 0.5*dy(j)
!            facexData(VELC_FACE_VAR,i,j,k) = -EXP(-2.0*dr_simTime)*COS(xedge)*SIN(ycell)
          enddo
        enddo
      enddo

      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, FACEY)
      do k=1, blkLimitsGC(HIGH,KAXIS)
        do j=1, blkLimitsGC(HIGH,JAXIS)
          do i=1, blkLimitsGC(HIGH,IAXIS)
!            xcell = coord(IAXIS) - bsize(IAXIS)/2.0 + real(i-NGUARD-1)*dx(i) + 0.5*dx(i)
!            yedge = coord(JAXIS) - bsize(JAXIS)/2.0 + real(j-NGUARD-1)*dy(j)
!            faceyData(VELC_FACE_VAR,i,j,k) = EXP(-2.0*dr_simTime)*COS(yedge)*SIN(xcell)
          enddo
        enddo
      enddo
            
      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, CENTER)
      do k=1, blkLimitsGC(HIGH,KAXIS)
        do j=1, blkLimitsGC(HIGH,JAXIS)
          do i=1, blkLimitsGC(HIGH,IAXIS)
!            xcell = coord(IAXIS) - bsize(IAXIS)/2.0 + real(i-NGUARD-1)*dx(i) + 0.5*dx(i)
!            ycell = coord(JAXIS) - bsize(JAXIS)/2.0 + real(j-NGUARD-1)*dy(j) + 0.5*dy(j)
!            solnData(TEMP_VAR,i,j,k) = -EXP(-2.0*dr_simTime)*COS(xcell)*COS(ycell)/2.0 + 0.5
          enddo
        enddo
      enddo

      solnData(PRES_VAR,:,:,:) = 0.0
      solnData(DELP_VAR,:,:,:) = 0.0
      solnData(DUST_VAR,:,:,:) = 0.0
      solnData(TVIS_VAR,:,:,:) = 0.0

      facexData(RHDS_FACE_VAR,:,:,:) = 0.0
      faceyData(RHDS_FACE_VAR,:,:,:) = 0.0
    
#if NDIM == 3
      facezData(VELC_FACE_VAR,:,:,:) = 0.0
      facezData(RHDS_FACE_VAR,:,:,:) = 0.0
#endif      




    ! taylor-green flow (Modified Temp)
    case (21)
      
      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, CENTER)
      do k=1, blkLimitsGC(HIGH,KAXIS)
        do j=1, blkLimitsGC(HIGH,JAXIS)
          do i=1, blkLimitsGC(HIGH,IAXIS)
             facexData(VELC_FACE_VAR,i,j,k) = -EXP(-2.0*dr_simTime)*COS(xedge(i))*SIN(ycell(j))
          enddo
        enddo
      enddo

      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, CENTER)
      do k=1, blkLimitsGC(HIGH,KAXIS)
        do j=1, blkLimitsGC(HIGH,JAXIS)
          do i=1, blkLimitsGC(HIGH,IAXIS)
            faceyData(VELC_FACE_VAR,i,j,k) = EXP(-2.0*dr_simTime)*COS(yedge(j))*SIN(xcell(i))
          enddo
        enddo
      enddo
            
      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, CENTER)
      do k=1, blkLimitsGC(HIGH,KAXIS)
        do j=1, blkLimitsGC(HIGH,JAXIS)
          do i=1, blkLimitsGC(HIGH,IAXIS)
            solnData(TEMP_VAR,i,j,k) = -EXP(-2.0*dr_simTime)*SIN(xcell(i)/2.)*SIN(ycell(j)/2.)/2.0 + 0.5
          enddo
        enddo
      enddo

      solnData(PRES_VAR,:,:,:) = 0.0
      solnData(DELP_VAR,:,:,:) = 0.0
      solnData(DUST_VAR,:,:,:) = 0.0
      solnData(TVIS_VAR,:,:,:) = 0.0

      facexData(RHDS_FACE_VAR,:,:,:) = 0.0
      faceyData(RHDS_FACE_VAR,:,:,:) = 0.0

#if NDIM == 3
      facezData(VELC_FACE_VAR,:,:,:) = 0.0
      facezData(RHDS_FACE_VAR,:,:,:) = 0.0
#endif      




    
    ! null initialization  
    case default

      solnData(PRES_VAR,:,:,:) = 0.0
      solnData(DELP_VAR,:,:,:) = 0.0
      solnData(DUST_VAR,:,:,:) = 0.0
      solnData(TVIS_VAR,:,:,:) = 0.0
      solnData(TEMP_VAR,:,:,:) = 0.0

      facexData(VELC_FACE_VAR,:,:,:) = 0.0
      faceyData(VELC_FACE_VAR,:,:,:) = 0.0
      facexData(RHDS_FACE_VAR,:,:,:) = 0.0
      faceyData(RHDS_FACE_VAR,:,:,:) = 0.0

#if NDIM == 3
      facezData(VELC_FACE_VAR,:,:,:) = 0.0
      facezData(RHDS_FACE_VAR,:,:,:) = 0.0
#endif      




  
  end select

  call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  call Grid_releaseBlkPtr(blockID,facexData,FACEX)
  call Grid_releaseBlkPtr(blockID,faceyData,FACEY)

#if NDIM == 3
  call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#endif


  return

111    format (i4,3x,i4)
112    format (3(3x,e12.4))

end subroutine Simulation_initBlock
