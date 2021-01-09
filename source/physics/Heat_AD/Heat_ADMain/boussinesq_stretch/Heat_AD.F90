subroutine Heat_AD(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)

  use Heat_AD_data, only: ht_gama_coeff, ht_rho_coeff, &
                          ht_gama, ht_rho
  use Heat_AD_interface, only: Heat_Solve, ht_rhs2d, ht_rhs3d
  use Grid_interface,    only: Grid_getCellMetrics,    &
                               Grid_getBlkIndexLimits, &
                               Grid_getBlkPtr,         &
                               Grid_releaseBlkPtr,     &
                               Grid_fillGuardCells

  implicit none

#include "constants.h"
#include "Heat_AD.h"
#include "Flash.h"   

  include "Flash_mpi.h"

  integer, INTENT(INOUT) :: blockCount
  integer, INTENT(IN) :: sweepOrder
  integer, INTENT(INOUT) :: blockList(MAXBLOCKS)
  real,    INTENT(IN) :: timeEndAdv, dt, dtOld

  integer ::  blockID,lb,ist,itmx
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)
  real, pointer, dimension(:,:,:,:) :: solnData, facexData, faceyData, facezData
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: Trhs
  real, dimension(GRID_IHI_GC,3,blockCount) :: iMetrics
  real, dimension(GRID_JHI_GC,3,blockCount) :: jMetrics
  real, dimension(GRID_KHI_GC,3,blockCount) :: kMetrics

  Trhs = 0.

  do lb = 1,blockCount
    blockID = blockList(lb)

    call Grid_getCellMetrics(IAXIS,blockID,LEFT_EDGE, .true.,iMetrics(:,LEFT_EDGE,lb), GRID_IHI_GC)
    call Grid_getCellMetrics(IAXIS,blockID,CENTER,    .true.,iMetrics(:,CENTER,lb),    GRID_IHI_GC)
    call Grid_getCellMetrics(IAXIS,blockID,RIGHT_EDGE,.true.,iMetrics(:,RIGHT_EDGE,lb),GRID_IHI_GC)

    call Grid_getCellMetrics(JAXIS,blockID,LEFT_EDGE, .true.,jMetrics(:,LEFT_EDGE,lb), GRID_JHI_GC)
    call Grid_getCellMetrics(JAXIS,blockID,CENTER,    .true.,jMetrics(:,CENTER,lb),    GRID_JHI_GC)
    call Grid_getCellMetrics(JAXIS,blockID,RIGHT_EDGE,.true.,jMetrics(:,RIGHT_EDGE,lb),GRID_JHI_GC)

    call Grid_getCellMetrics(KAXIS,blockID,LEFT_EDGE, .true.,kMetrics(:,LEFT_EDGE,lb), GRID_KHI_GC)
    call Grid_getCellMetrics(KAXIS,blockID,CENTER,    .true.,kMetrics(:,CENTER,lb),    GRID_KHI_GC)
    call Grid_getCellMetrics(KAXIS,blockID,RIGHT_EDGE,.true.,kMetrics(:,RIGHT_EDGE,lb),GRID_KHI_GC)
  end do

  ! 3rd order Runge Kutta Coefficents
  if ( .true. ) then

    ht_gama_coeff(1) = 8./15.
    ht_gama_coeff(2) = 5./12.
    ht_gama_coeff(3) = 3./4.

    ht_rho_coeff(1) = 0.
    ht_rho_coeff(2) = -17./60.
    ht_rho_coeff(3) = -5./12.

    itmx = 3

  ! 2nd order Adams Bashforth Coefficents (const dt)
  elseif ( .false. ) then

    ht_gama_coeff(1) = 1.5
    ht_gama_coeff(2) = 0.0
    ht_gama_coeff(3) = 0.0

    ht_rho_coeff(1) = -0.5
    ht_rho_coeff(2) = 0.0
    ht_rho_coeff(3) = 0.0

    itmx = 1

  ! 2nd order Adams Bashforth Coefficents (variable dt)
  else

    ht_gama_coeff(1) = 1.5
    ht_gama_coeff(2) = 0.0
    ht_gama_coeff(3) = 0.0

    ht_rho_coeff(1) = -0.5
    ht_rho_coeff(2) = 0.0
    ht_rho_coeff(3) = 0.0

    itmx = 1

  endif

  ! Time integration loop
  do ist = 1, itmx

    ht_gama = ht_gama_coeff(ist)
    ht_rho = ht_rho_coeff(ist)

    do lb = 1,blockCount
      blockID = blockList(lb)

      ! Get Blocks internal limits indexes:
      call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

      ! Point to blocks center and face vars:
      call Grid_getBlkPtr(blockID,solnData,CENTER)
      call Grid_getBlkPtr(blockID,facexData,FACEX)
      call Grid_getBlkPtr(blockID,faceyData,FACEY)

#if NDIM == 3
      call Grid_getBlkPtr(blockID,facezData,FACEZ)

      call ht_rhs3d(solnData(TEMP_VAR,:,:,:), Trhs,    &
                    facexData(VOLD_FACE_VAR,:,:,:),    &
                    faceyData(VOLD_FACE_VAR,:,:,:),    &
                    facezData(VOLD_FACE_VAR,:,:,:),    &
                    iMetrics(:,:,lb),                  &
                    jMetrics(:,:,lb),                  &
                    kMetrics(:,:,lb),                  &
           blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS), &
           blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS), &
           blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS))

#elif NDIM == 2

      call ht_rhs2d(solnData(TEMP_VAR,:,:,:), Trhs,    &
                    facexData(VOLD_FACE_VAR,:,:,:),    &
                    faceyData(VOLD_FACE_VAR,:,:,:),    &
                    iMetrics(:,:,lb),jMetrics(:,:,lb), &
           blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS), &
           blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS))

#endif

      call Heat_Solve(solnData(TEMP_VAR,:,:,:), Trhs,  &
                      solnData(TOLD_VAR,:,:,:), dt,    &
           blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS), &
           blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS), &
           blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS), &
                      ht_gama, ht_rho)

      solnData(TOLD_VAR,:,:,:) = Trhs(:,:,:)

      call Grid_releaseBlkPtr(blockID,solnData,CENTER)
      call Grid_releaseBlkPtr(blockID,facexData,FACEX)
      call Grid_releaseBlkPtr(blockID,faceyData,FACEY)

#if NDIM == 3
      call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

    end do
  
  gcMask = .FALSE.
  gcMask(TEMP_VAR) = .true.
  call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)

  end do

end subroutine Heat_AD
