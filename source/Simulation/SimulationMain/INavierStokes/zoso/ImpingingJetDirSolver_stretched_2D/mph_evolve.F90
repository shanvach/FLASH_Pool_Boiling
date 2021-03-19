subroutine mph_evolve(blockCount, blockList, timeEndAdv,dt,dtOld,sweepOrder,mph_flag) 

  ! Following routine is written by Akash
  ! Actual calls written by Shizao and Keegan
  ! This subroutine decouples Multiphase calls from ins_ab2rk3_VD 

!#define NUCLEATE_BOILING
#include "Flash.h"

  ! Modules Use:
#ifdef FLASH_GRID_PARAMESH
  use physicaldata, ONLY : interp_mask_unk_res,      &
                           interp_mask_facex_res,    &
                           interp_mask_facey_res,    &
                           interp_mask_facez_res,    &
                           interp_mask_unk,      &
                           interp_mask_facex,    &
                           interp_mask_facey,    &
                           interp_mask_facez
  use workspace, ONLY :    interp_mask_work
#endif    

  use Grid_interface, ONLY : Grid_getDeltas,         &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getCellMetrics,    &
                             Grid_getCellCoords,     &
                             Grid_getBlkIndexLimits, &
                             Grid_fillGuardCells,    &
                             Grid_getBlkBoundBox,Grid_getBlkCenterCoords

  use Multiphase_data, only: mph_rho1,mph_rho2,mph_sten,mph_crmx,mph_crmn, &
                             mph_vis1,mph_vis2,mph_lsit, mph_inls,mph_meshMe, &
                             mph_jet_vel,mph_jet_src,mph_jet_tstamp,mph_jet_twait, &
                             mph_jet_period,mph_dist_flag,mph_vel_scalar

  use mph_interface, only : mph_KPDcurvature2DAB, mph_KPDcurvature2DC, &
                            mph_KPDadvectWENO3, mph_KPDlsRedistance,  &
                            mph_KPDcurvature3DAB, mph_KPDcurvature3DC,&
                            mph_KPDadvectWENO3_3D,mph_KPDlsRedistance_3D

  use Timers_interface, ONLY : Timers_start, Timers_stop

  use Driver_data, ONLY : dr_nstep,dr_simTime

  use ins_interface, only: ins_fluxfixRho1,ins_fluxfixRho2

  use IncompNS_data, only: ins_gravY

  use Simulation_data, only: sim_jet_depth, sim_jet_x, &
                             sim_jet_z, sim_free_surface, sim_yMin,&
                             sim_xMax

  implicit none

#include "constants.h"
#include "IncompNS.h"
  include "Flash_mpi.h"

  ! Arugment List
  integer, intent(in) :: sweepOrder
  integer, INTENT(INOUT) :: blockCount
  integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList
  real,    INTENT(IN) :: timeEndAdv,dt,dtOld
  integer, intent(in) :: mph_flag

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  real, dimension(2,MDIM) :: boundBox

  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

  integer :: lb,blockID,ii,i,j,k,ierr

  real, dimension(GRID_IHI_GC,3,blockCount) :: dx
  real, dimension(GRID_JHI_GC,3,blockCount) :: dy
  real, dimension(GRID_KHI_GC,3,blockCount) :: dz

  real bsize(MDIM),coord(MDIM),xCoord(GRID_IHI_GC),xFace(GRID_IHI_GC+1)

  real del(MDIM)

  !kpd
  real :: lsDT,lsT,minCellDiag
  real :: volSum,volSumAll,eps
  
  !- kpd - For Overall Solver Timer... 
  real :: t_startMP1,t_stopMP1,t_startMP2,t_startMP2a,t_stopMP2

  integer :: t

  integer :: listofBlocks(MAXBLOCKS)
  integer :: count
  integer :: intval,nxc,nyc,nzc

  real :: xcell,ycell,zcell,theta,interface_loc

  real :: ycell_plus, dfun_y, dfun_y_plus, dxc, dlen, sgnr, sgnl, jet

  real :: distLeft, distRight, dist, xLeft, xRight

  real :: free_surface_local
  
  real :: sgn = 1.0

   do lb = 1,blockCount

     blockID = blockList(lb)

  ! Get blk cell metrics by direction from Grid Unit 
    
     call Grid_getCellMetrics(IAXIS,blockID,LEFT_EDGE, .true.,dx(:,LEFT_EDGE,lb), GRID_IHI_GC) 
     call Grid_getCellMetrics(IAXIS,blockID,CENTER,    .true.,dx(:,CENTER,lb),    GRID_IHI_GC) 
     call Grid_getCellMetrics(IAXIS,blockID,RIGHT_EDGE,.true.,dx(:,RIGHT_EDGE,lb),GRID_IHI_GC) 

     call Grid_getCellMetrics(JAXIS,blockID,LEFT_EDGE, .true.,dy(:,LEFT_EDGE,lb), GRID_JHI_GC) 
     call Grid_getCellMetrics(JAXIS,blockID,CENTER,    .true.,dy(:,CENTER,lb),    GRID_JHI_GC) 
     call Grid_getCellMetrics(JAXIS,blockID,RIGHT_EDGE,.true.,dy(:,RIGHT_EDGE,lb),GRID_JHI_GC) 


     call Grid_getCellMetrics(KAXIS,blockID,LEFT_EDGE, .true.,dz(:,LEFT_EDGE,lb), GRID_KHI_GC) 
     call Grid_getCellMetrics(KAXIS,blockID,CENTER,    .true.,dz(:,CENTER,lb),    GRID_KHI_GC) 
     call Grid_getCellMetrics(KAXIS,blockID,RIGHT_EDGE,.true.,dz(:,RIGHT_EDGE,lb),GRID_KHI_GC) 

  enddo  

if(mph_flag == 1) then

!kpd - Level Set Initialization...
  if (dr_nstep .eq. 1) then

    !#############################################################################
    !-----------------------------------------------------------------------------
    !-kpd - Fill Guardcells for the distance function before curvature is
    !computed 
    !       This is done for first iteration only (filled at end of time step) 
    !-----------------------------------------------------------------------------
    
    gcMask = .FALSE.
    gcMask(DFUN_VAR) = .TRUE.
#ifdef FLASH_GRID_PARAMESH
    intval = 1
    !intval = 2
    interp_mask_unk = intval;   interp_mask_unk_res = intval;
    interp_mask_work= intval;
#endif

    call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask,selectBlockType=ACTIVE_BLKS)

   !*********************************************************************************************************
   !- kpd - Level Set Distance Function Initialization (if needed)
   !******************************************
   !*********************************************************************************************************
   lsT  = 0.0
   do ii = 1,mph_inls

     !------------------------------
     !- kpd - Level set redistancing 
     !------------------------------

     t = dt

     do lb = 1,blockCount
        blockID = blockList(lb)

        ! Get blocks dx, dy ,dz:
        !call Grid_getDeltas(blockID,del)

        ! Get Blocks internal limits indexes:
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

        ! Point to blocks center and face vars:
        call Grid_getBlkPtr(blockID,solnData,CENTER)
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
        call Grid_getBlkPtr(blockID,facezData,FACEZ)

        !--------------------------------------------
        ! Call DFUN re-initialization routine for 3D:
        !--------------------------------------------
        lsDT = MIN(10.0*dt,0.001)
        !minCellDiag = SQRT(del(DIR_X)**2.+del(DIR_Y)**2.+del(DIR_Z)**2.)
        minCellDiag = 0.001 !SQRT((SQRT(del(DIR_X)**2.+del(DIR_Y)**2.))**2.+del(DIR_Z)**2.)
        if ( ii .eq. mph_inls .AND. lb .eq. 1 .AND. mph_meshMe .eq. 0) then
           print*,"Level Set Initialization Iteration # ",ii,minCellDiag,lsDT
        end if

        if (ii.eq.1) solnData(AAJUNK_VAR,:,:,:) = solnData(DFUN_VAR,:,:,:)

        call mph_KPDlsRedistance_3D(solnData(DFUN_VAR,:,:,:), &
                          facexData(VELC_FACE_VAR,:,:,:), &
                          faceyData(VELC_FACE_VAR,:,:,:), &
                          facezData(VELC_FACE_VAR,:,:,:), &
                          dx(:,:,blockID),dy(:,:,blockID),dz(:,:,blockID),  &
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS), &
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS), &
                          blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS), &
                          solnData(AAJUNK_VAR,:,:,:), lsDT, minCellDiag )

#elif NDIM == 2


        !--------------------------------------------
        ! Call DFUN re-initialization routine for 2D:
        !--------------------------------------------
        !lsDT = MIN(10.0*dt,0.001)
        minCellDiag = 0.001 !SQRT(del(DIR_X)**2.+del(DIR_Y)**2.)

        lsDT = minCellDiag/2.0d0
        if ( ii .eq. mph_inls .AND. lb .eq. 1 .AND. mph_meshMe .eq. 0) then
           print*,"Level Set Initialization Iteration # ",ii,minCellDiag,lsDT
        end if

        if (ii.eq.1) solnData(AAJUNK_VAR,:,:,:) = solnData(DFUN_VAR,:,:,:)

        call mph_KPDlsRedistance(solnData(DFUN_VAR,:,:,:), &
                          facexData(VELC_FACE_VAR,:,:,:),  &
                          faceyData(VELC_FACE_VAR,:,:,:),  &
                          dx(:,:,blockID),dy(:,:,blockID),  &
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS), &
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS), &
                          solnData(AAJUNK_VAR,:,:,:), lsDT, blockID,minCellDiag)

#endif

        ! Release pointers:
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM ==3
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
     enddo    !do lb = 1,blockCount

    !*********************************************************************************************************
    !-kpd - Fill distance function guard cells after each re-initialization to
    !communicate updates
    gcMask = .FALSE.
    gcMask(DFUN_VAR) = .TRUE.
#ifdef FLASH_GRID_PARAMESH
    intval = 1
    !intval = 2
    interp_mask_unk = intval;   interp_mask_unk_res = intval;
    interp_mask_work= intval;
#endif
    call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)
    !*********************************************************************************************************
    !call mph_bcLevelSet(0.0d0,0.0d0)

      lsT = lsT + lsDT

     if(mod(ii,1000) == 0) then
       write(*,*) 'Output DFUN:', mph_meshMe, ii, lsT, lsDT
       call outtotecplot(mph_meshMe,lsT,lsDT,ii,ii, &
                         0.0,blockList,blockCount,1)
     endif

   end do  ! End do: ii=1,inls
           !==================

   end if  ! End if: dr_nstep = 1
           !=====================

!************************************************************************************************************
!************************************************************************************************************

    !-----------------------------------------------------
    !- kpd - Loop through current block for curvature 2dA
    !-----------------------------------------------------

    do lb = 1,blockCount
     blockID = blockList(lb)
!    do lb = 1,count
!     blockID = listofBlocks(lb)

     !----------------------------------------------------------
     !- kpd - Get Block Information...
     !----------------------------------------------------------
     call Grid_getBlkBoundBox(blockId,boundBox)
     call Grid_getBlkBoundBox(blockId,boundBox)
     bsize(:) = boundBox(2,:) - boundBox(1,:)
     call Grid_getBlkCenterCoords(blockId,coord)

     ! Get block's dx, dy ,dz:
     !call Grid_getDeltas(blockID,del)

     ! Get block's internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     ! Point to block's center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     !----------------------------------------------------------

     !----------------------------------------------------------
     !!- kpd - Screen output for AMR testing
     !if (ins_nstep .eq. 1)print*,"KPD ins_ab2rk3
     !bID=",blockID,"blockCount=",blockCount,coord(1),coord(2),coord(3)
     !if (ins_nstep .eq. 1)print*,"KPD ins_ab2rk3
     !bID=",blockID,coord(1),coord(2),coord(3)
     !----------------------------------------------------------
#if NDIM == 2

     !----------------------------------------------------------
     !- kpd - Call 2-D curvature Routine:
     !----------------------------------------------------------
     ! Akash - Modified call to compute specific heat and thermal conductivity

     call mph_KPDcurvature2DAB(solnData(DFUN_VAR,:,:,:),               &
                           solnData(CURV_VAR,:,:,:),                   &
                           facexData(RH1F_FACE_VAR,:,:,:),             &
                           facexData(RH2F_FACE_VAR,:,:,:),             &
                           faceyData(RH1F_FACE_VAR,:,:,:),             &
                           faceyData(RH2F_FACE_VAR,:,:,:),             &
                           solnData(PFUN_VAR,:,:,:),                   &
                           solnData(SIGP_VAR,:,:,:),                   &
                           facexData(SIGM_FACE_VAR,:,:,:),             &
                           faceyData(SIGM_FACE_VAR,:,:,:),             &
                           dx(:,:,blockID),dy(:,:,blockID),mph_rho1,mph_rho2,    &
                           mph_sten,mph_crmx,mph_crmn,                 &
                           blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS), &
                           blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS), &
                           solnData(VISC_VAR,:,:,:),mph_vis1,mph_vis2)
     !------------------------------------------------------------------

#elif NDIM ==3

        call Grid_getBlkPtr(blockID,facezData,FACEZ)

        !----------------------------------------------------------
        !- kpd - Call curvature3DAB Routine to compute curvature,
        !           phase function, and face densities:
        !----------------------------------------------------------
        call mph_KPDcurvature3DAB(solnData(DFUN_VAR,:,:,:),            &
                           solnData(CURV_VAR,:,:,:),                   &
                           dx(:,:,blockID),dy(:,:,blockID), dz(:,:,blockID),&
                           blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS), &
                           blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS), &
                           blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS), &
                           facexData(RH1F_FACE_VAR,:,:,:),             &
                           facexData(RH2F_FACE_VAR,:,:,:),             &
                           faceyData(RH1F_FACE_VAR,:,:,:),             &
                           faceyData(RH2F_FACE_VAR,:,:,:),             &
                           facezData(RH1F_FACE_VAR,:,:,:),             &
                           facezData(RH2F_FACE_VAR,:,:,:),             &
                           solnData(PFUN_VAR,:,:,:),                   &
                           mph_rho1,mph_rho2,                          &
                           solnData(VISC_VAR,:,:,:),mph_vis1,mph_vis2)
        !----------------------------------------------------------

#endif

     !-----------------------------------------------
     ! Release pointers:
     !-----------------------------------------------
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM ==3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
     !-----------------------------------------------

    enddo

  !--------------------------------------------------------------------------
  !- kpd - Implemented for variable density. Fill multiphase Guard Cell data.
  ! -------------------------------------------------------------------------
  gcMask = .FALSE.
  gcMask(PFUN_VAR) = .TRUE.                                ! Phase Function
  gcMask(CURV_VAR) = .TRUE.                                ! Curvature
  gcMask(VISC_VAR) = .TRUE.                                ! Viscosity

  gcMask(NUNK_VARS+RH1F_FACE_VAR) = .TRUE.                 ! rho1x
  gcMask(NUNK_VARS+1*NFACE_VARS+RH1F_FACE_VAR) = .TRUE.    ! rho1y
  gcMask(NUNK_VARS+RH2F_FACE_VAR) = .TRUE.                 ! rho2x
  gcMask(NUNK_VARS+1*NFACE_VARS+RH2F_FACE_VAR) = .TRUE.    ! rho2y

#if NDIM == 3
  gcMask(NUNK_VARS+2*NFACE_VARS+RH1F_FACE_VAR) = .TRUE.    ! rho1z
  gcMask(NUNK_VARS+2*NFACE_VARS+RH2F_FACE_VAR) = .TRUE.    ! rho2z
#endif

#ifdef FLASH_GRID_PARAMESH
  intval = 1
  !intval = 2
  interp_mask_unk = intval;   interp_mask_unk_res = intval;
  interp_mask_work= intval;
  interp_mask_facex = intval; interp_mask_facex_res = intval;
  interp_mask_facey = intval; interp_mask_facey_res = intval;
  interp_mask_facez = intval; interp_mask_facez_res = intval;
#endif
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask,selectBlockType=ACTIVE_BLKS)
       !maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)

!***********************************************************************************************
!***********************************************************************************************
!***********************************************************************************************

else if(mph_flag == 0) then

    !-----------------------------------------------------
    !- kpd - Loop through current block for curvature 2dC
    !-----------------------------------------------------

    do lb = 1,blockCount
     blockID = blockList(lb)
!    do lb = 1,count
!     blockID = listofBlocks(lb)

     ! Get blocks dx, dy ,dz:
     !call Grid_getDeltas(blockID,del)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

#if NDIM == 2

     !----------------------------------------------------------
     !- kpd - Call 2-D curvature Routine:
     !----------------------------------------------------------
     call mph_KPDcurvature2DC(solnData(DFUN_VAR,:,:,:), &
                          solnData(CURV_VAR,:,:,:), &
                          facexData(RH1F_FACE_VAR,:,:,:), &
                          facexData(RH2F_FACE_VAR,:,:,:), &
                          faceyData(RH1F_FACE_VAR,:,:,:), &
                          faceyData(RH2F_FACE_VAR,:,:,:), &
                          solnData(PFUN_VAR,:,:,:), &
                          solnData(SIGP_VAR,:,:,:), &
                          solnData(SPCD_VAR,:,:,:),&
                          facexData(SIGM_FACE_VAR,:,:,:), &
                          faceyData(SIGM_FACE_VAR,:,:,:), &
                          facexData(SIGC_FACE_VAR,:,:,:), &
                          faceyData(SIGC_FACE_VAR,:,:,:), & 
                          dx(:,:,blockID),dy(:,:,blockID),mph_rho1,mph_rho2, &
                          mph_sten,mph_crmx,mph_crmn, &
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),blockID)

#elif NDIM == 3 
        call Grid_getBlkPtr(blockID,facezData,FACEZ)

        !----------------------------------------------------------
        !- kpd - Call curvature3DC Routine to compute interfacial
        !           densities, and momentum/poisson jumps:
        !----------------------------------------------------------

        call mph_KPDcurvature3DC( solnData(DFUN_VAR,:,:,:)  , &
                           solnData(CURV_VAR,:,:,:)         , &
                           facexData(RH1F_FACE_VAR,:,:,:)   , &
                           facexData(RH2F_FACE_VAR,:,:,:)   , &
                           faceyData(RH1F_FACE_VAR,:,:,:)   , &
                           faceyData(RH2F_FACE_VAR,:,:,:)   , &
                           solnData(PFUN_VAR,:,:,:)         , &
                           solnData(SIGP_VAR,:,:,:)         , &
                           solnData(SPCD_VAR,:,:,:)         , &
                           facexData(SIGM_FACE_VAR,:,:,:)   , &
                           faceyData(SIGM_FACE_VAR,:,:,:)   , &
                           dx(:,:,blockID),dy(:,:,blockID)  , &
                           mph_rho1,mph_rho2                , &
                           mph_sten                         , &
                           blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                           blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                           dz(:,:,blockID)                  , &
                           blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
                           facezData(RH1F_FACE_VAR,:,:,:)   , &
                           facezData(RH2F_FACE_VAR,:,:,:)   , &
                           facezData(SIGM_FACE_VAR,:,:,:))

#endif
     !-----------------------------------------------
     ! Release pointers:
     !-----------------------------------------------
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM ==3
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
     !-----------------------------------------------

    enddo

  !--------------------------------------------------------------------------
  !- kpd - Implemented for variable density. Fill multiphase Guard Cell data.
  ! -------------------------------------------------------------------------
  gcMask = .FALSE.
  gcMask(NUNK_VARS+RH1F_FACE_VAR) = .TRUE.                 ! rho1x
  gcMask(NUNK_VARS+1*NFACE_VARS+RH1F_FACE_VAR) = .TRUE.    ! rho1y
  gcMask(NUNK_VARS+RH2F_FACE_VAR) = .TRUE.                 ! rho2x
  gcMask(NUNK_VARS+1*NFACE_VARS+RH2F_FACE_VAR) = .TRUE.    ! rho2y

  gcMask(SIGP_VAR) = .TRUE.                                ! Poisson Jump
  gcMask(NUNK_VARS+SIGM_FACE_VAR) = .TRUE.                 ! Momentum Jump X
  gcMask(NUNK_VARS+1*NFACE_VARS+SIGM_FACE_VAR) = .TRUE.    ! Momentum Jump Y

  gcMask(SPCD_VAR) = .TRUE.                                ! Poisson Jump 
  gcMask(NUNK_VARS+SIGC_FACE_VAR) = .TRUE.                 ! Momentum Jump X
  gcMask(NUNK_VARS+1*NFACE_VARS+SIGC_FACE_VAR) = .TRUE.    ! Momentum Jump Y


  gcMask(CURV_VAR) = .TRUE.

#if NDIM == 3
  gcMask(NUNK_VARS+2*NFACE_VARS+RH1F_FACE_VAR) = .TRUE.    ! rho1z
  gcMask(NUNK_VARS+2*NFACE_VARS+RH2F_FACE_VAR) = .TRUE.    ! rho2z
  gcMask(NUNK_VARS+2*NFACE_VARS+SIGM_FACE_VAR) = .TRUE.    ! Momentum Jump Z
#endif
#ifdef FLASH_GRID_PARAMESH
  intval = 1
  !intval = 2
  interp_mask_unk = intval;   interp_mask_unk_res = intval;
  interp_mask_work= intval;
  interp_mask_facex = intval; interp_mask_facex_res = intval;
  interp_mask_facey = intval; interp_mask_facey_res = intval;
  interp_mask_facez = intval; interp_mask_facez_res = intval;
#endif

  !print*,"KPD - Filling Density Guard Cells."
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask,selectBlockType=ACTIVE_BLKS)

  nxc = NXB + NGUARD + 1
  nyc = NYB + NGUARD + 1
  nzc = NZB + NGUARD + 1  

#ifdef FLASH_GRID_PARAMESH
  call ins_fluxfixRho1(NGUARD,nxc,nyc,nzc,nxc-1,nyc-1,nzc-1,&
                   blockCount,blockList)
#endif

#ifdef FLASH_GRID_PARAMESH
  call ins_fluxfixRho2(NGUARD,nxc,nyc,nzc,nxc-1,nyc-1,nzc-1,&
                   blockCount,blockList)
#endif

  jet = 1.0

  if(mph_dist_flag .and. dr_simTime .lt. (mph_jet_tstamp + mph_jet_period)) then
 
    mph_vel_scalar = 1.0 !+ 2.5*abs(sin((dr_simTime - mph_jet_tstamp)*3.14/mph_jet_period))

  else
   
    if(mph_dist_flag) mph_jet_tstamp = dr_simTime

    mph_dist_flag = .false.
    mph_vel_scalar = 1.0

    if(dr_simTime .gt. (mph_jet_tstamp + mph_jet_twait)) then
        mph_dist_flag = .true.
        mph_jet_tstamp = dr_simTime
    end if
    
  end if

  if(mph_meshMe .eq. MASTER_PE) print *,"Jet Velocity:", dr_simTime,  mph_vel_scalar

  mph_jet_vel(:,:,:) = 0.0

  do lb = 1,blockCount

     blockID = blockList(lb)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     call Grid_getBlkBoundBox(blockId,boundBox)
     bsize(:) = boundBox(2,:) - boundBox(1,:)

     !call Grid_getBlkCenterCoords(blockId,coord)

     !call Grid_getDeltas(blockID,del)

     call Grid_getBlkPtr(blockID,solnData,CENTER)

     call Grid_getCellCoords(IAXIS,lb,CENTER,.true.,xCoord,blkLimitsGC(HIGH,IAXIS)-1)
     call Grid_getCellCoords(IAXIS,lb,FACES,.true.,xFace,blkLimitsGC(HIGH,IAXIS))

     k = 1
   
     sim_jet_x = 0.0

     sim_jet_z = 0.0

     do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
#if NDIM == 3
     do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
#endif

       !xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
       !        real(i - NGUARD - 1)/dx(i,CENTER,blockID) +   &
       !        0.5/dx(i,CENTER,blockID)
       xcell = xCoord(i)
       xLeft = xFace(i)
       xRight = xFace(i+1)
       zcell = 0.0

#if NDIM == 3
       zcell = coord(KAXIS) - bsize(KAXIS)/2.0 +   &
               real(k - NGUARD - 1)*(1/dz(k,CENTER,blockID)) +   &
               0.5*(1/dz(k,CENTER,blockID))
#endif

       dxc = 1/dx(i,CENTER,blockID)
       dist = sqrt((xcell-sim_jet_x)**2+(zcell-sim_jet_z)**2)
       distLeft = sqrt((xLeft-sim_jet_x)**2+(zcell-sim_jet_z)**2)
       distRight = sqrt((xRight-sim_jet_x)**2+(zcell-sim_jet_z)**2)

       sgnr = sign(sgn,0.5-distRight)
       sgnl = sign(sgn,0.5-distLeft)

       dlen = ((sgnr-sgnl)*0.5)*sgnr*(((0.5*jet-distRight)/dxc)*((sgnr+1)*0.5)+&
                                      ((0.5*jet-distLeft )/dxc)*((sgnl+1)*0.5))+&
              0.25*(sgnr*sgnl+sgnl)*(sgnr*sgnl+sgnr)
              

       mph_jet_vel(i,k,blockID) = mph_vel_scalar*dlen
       !print*,mph_jet_vel(i,k,blockID),sgnr,sgnl,dlen
       !if(distRight .lt. 0.5 .and. distLeft .gt. 0.5) then
       !   dlen = 0.5-distRight
       !   mph_jet_vel(i,k,blockID) = mph_vel_scalar*(dlen/dxc)
       !elseif(distRight .gt. 0.5 .and. distLeft .lt. 0.5) then
       !   dlen = 0.5-distLeft
       !   mph_jet_vel(i,k,blockID) = mph_vel_scalar*(dlen/dxc)
       !elseif(distRight .lt. 0.5 .and. distLeft .lt. 0.5) then
       !   dlen = dxc
       !   mph_jet_vel(i,k,blockID) = mph_vel_scalar*(dlen/dxc)
       !else
       !   mph_jet_vel(i,k,blockID) = 0.0
       !endif
 
       mph_jet_src(i,k,blockID)  = sqrt((xcell-sim_jet_x)**2+(zcell-sim_jet_z)**2) - 0.5

#if NDIM == 3
     end do
#endif

     end do
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  end do

end if

end subroutine
