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
                             Grid_getBlkIndexLimits, &
                             Grid_fillGuardCells,    &
                             Grid_getBlkBoundBox,Grid_getBlkCenterCoords

  use Multiphase_data, only: mph_rho1,mph_rho2,mph_sten,mph_crmx,mph_crmn, &
                             mph_vis1,mph_vis2,mph_lsit, mph_inls,mph_meshMe,&
                             mph_thco1,mph_thco2,mph_cp1,mph_cp2 ! Akash

  use mph_interface, only : mph_KPDcurvature2DAB, mph_KPDcurvature2DC, &
                            mph_KPDadvectWENO3, mph_KPDlsRedistance,  &
                            mph_KPDcurvature3DAB, mph_KPDcurvature3DC,&
                            mph_KPDadvectWENO3_3D,mph_KPDlsRedistance_3D,&
                            mph_getSmearedProperties2D,mph_getSmearedProperties3D

  use Timers_interface, ONLY : Timers_start, Timers_stop

  use Driver_data, ONLY : dr_nstep

  use ins_interface, only: ins_fluxfixRho1,ins_fluxfixRho2

  use Multiphase_data, only: mph_bcFlag

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

  integer :: lb,blockID,ii,i,j,k

  real bsize(MDIM),coord(MDIM)

  real del(MDIM)

  !kpd
  real :: lsDT,lsT,minCellDiag
  real :: volSum,volSumAll
  
  !- kpd - For Overall Solver Timer... 
  real :: t_startMP1,t_stopMP1,t_startMP2,t_startMP2a,t_stopMP2

  integer :: t

  integer :: listofBlocks(MAXBLOCKS)
  integer :: count
  integer :: intval,nxc,nyc,nzc

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
        call Grid_getDeltas(blockID,del)

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
        minCellDiag = SQRT((SQRT(del(DIR_X)**2.+del(DIR_Y)**2.))**2.+del(DIR_Z)**2.)
        if ( ii .eq. mph_inls .AND. lb .eq. 1 .AND. mph_meshMe .eq. 0) then
           print*,"Level Set Initialization Iteration # ",ii,minCellDiag,lsDT
        end if

        if (ii.eq.1) solnData(AAJUNK_VAR,:,:,:) = solnData(DFUN_VAR,:,:,:)

        call mph_KPDlsRedistance_3D(solnData(DFUN_VAR,:,:,:), &
                          facexData(VELI_FACE_VAR,:,:,:), &
                          faceyData(VELI_FACE_VAR,:,:,:), &
                          facezData(VELI_FACE_VAR,:,:,:), &
                          del(DIR_X),del(DIR_Y),del(DIR_Z),  &
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS), &
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS), &
                          blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS), &
                          solnData(AAJUNK_VAR,:,:,:), lsDT, minCellDiag )

#elif NDIM == 2


        !--------------------------------------------
        ! Call DFUN re-initialization routine for 2D:
        !--------------------------------------------
        !lsDT = MIN(10.0*dt,0.001)
        minCellDiag = SQRT(del(DIR_X)**2.+del(DIR_Y)**2.)

        lsDT = minCellDiag/2.0d0
        if ( ii .eq. mph_inls .AND. lb .eq. 1 .AND. mph_meshMe .eq. 0) then
           print*,"Level Set Initialization Iteration # ",ii,minCellDiag,lsDT
        end if

        if (ii.eq.1) solnData(AAJUNK_VAR,:,:,:) = solnData(DFUN_VAR,:,:,:)

        call mph_KPDlsRedistance(solnData(DFUN_VAR,:,:,:), &
                          facexData(VELI_FACE_VAR,:,:,:),  &
                          faceyData(VELI_FACE_VAR,:,:,:),  &
                          del(DIR_X),del(DIR_Y),  &
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

  do lb = 1,blockCount
     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,facezData,FACEZ)

     !----------------------------------------------------------
     !- ML - Call mph_getSmearedProperties2D
     !----------------------------------------------------------

#if NDIM == 2

     call mph_getSmearedProperties2D(solnData(DFUN_VAR,:,:,:), &
                          solnData(PFUN_VAR,:,:,:),&
                          del(DIR_X),del(DIR_Y),mph_rho1,mph_rho2, &
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                          solnData(NRMX_VAR,:,:,:),&
                          solnData(NRMY_VAR,:,:,:),&
                          solnData(SMHV_VAR,:,:,:),&
                          solnData(SMRH_VAR,:,:,:))
#endif

#if NDIM == 3

     call mph_getSmearedProperties3D(solnData(DFUN_VAR,:,:,:), &
                          solnData(PFUN_VAR,:,:,:),&
                          del(DIR_X),del(DIR_Y),del(DIR_Z),mph_rho1,mph_rho2, &
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                          blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
                          solnData(NRMX_VAR,:,:,:),&
                          solnData(NRMY_VAR,:,:,:),&
                          solnData(NRMZ_VAR,:,:,:),&
                          solnData(SMHV_VAR,:,:,:),&
                          solnData(SMRH_VAR,:,:,:))

#endif

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  end do

  gcMask = .FALSE.
  gcMask(SMHV_VAR) = .TRUE.
  gcMask(SMRH_VAR) = .TRUE.
  
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask,selectBlockType=ACTIVE_BLKS)

 
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
     call Grid_getDeltas(blockID,del)

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
                           del(DIR_X),del(DIR_Y),mph_rho1,mph_rho2,    &
                           mph_sten,mph_crmx,mph_crmn,                 &
                           blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS), &
                           blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS), &
                           solnData(VISC_VAR,:,:,:),mph_vis1,mph_vis2, &
                           solnData(ALPH_VAR,:,:,:),&
                           mph_thco1,mph_thco2,mph_cp1,mph_cp2,&
                           solnData(NRMX_VAR,:,:,:),solnData(NRMY_VAR,:,:,:),&
                           solnData(MFLG_VAR,:,:,:),solnData(SMHV_VAR,:,:,:),solnData(SMRH_VAR,:,:,:))!,blockID)
     !------------------------------------------------------------------

#elif NDIM ==3

        call Grid_getBlkPtr(blockID,facezData,FACEZ)

        !----------------------------------------------------------
        !- kpd - Call curvature3DAB Routine to compute curvature,
        !           phase function, and face densities:
        !----------------------------------------------------------
        call mph_KPDcurvature3DAB(solnData(DFUN_VAR,:,:,:),            &
                           solnData(CURV_VAR,:,:,:),                   &
                           del(DIR_X),del(DIR_Y), del(DIR_Z),          &
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
                           solnData(VISC_VAR,:,:,:),mph_vis1,mph_vis2 ,&
                           solnData(ALPH_VAR,:,:,:),&
                           mph_thco1,mph_thco2,mph_cp1,mph_cp2,&
                           solnData(NRMX_VAR,:,:,:),solnData(NRMY_VAR,:,:,:),&
                           solnData(NRMZ_VAR,:,:,:),solnData(MFLG_VAR,:,:,:),&
                           solnData(SMHV_VAR,:,:,:))
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

  gcMask(ALPH_VAR) = .TRUE.                 ! Thermal Diffusivity - Akash
  gcMask(NRMX_VAR) = .TRUE. 
  gcMask(NRMY_VAR) = .TRUE.
  gcMask(NRMZ_VAR) = .TRUE.
  gcMask(MFLG_VAR) = .TRUE.

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
     call Grid_getDeltas(blockID,del)

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
                          facexData(SIGM_FACE_VAR,:,:,:), &
                          faceyData(SIGM_FACE_VAR,:,:,:), &
                          del(DIR_X),del(DIR_Y),mph_rho1,mph_rho2, &
                          mph_sten,mph_crmx,mph_crmn, &
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                          mph_thco1,mph_thco2,mph_cp1,mph_cp2,solnData(MDOT_VAR,:,:,:))!,blockID)

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
                           facexData(SIGM_FACE_VAR,:,:,:)   , &
                           faceyData(SIGM_FACE_VAR,:,:,:)   , &
                           del(DIR_X),del(DIR_Y)            , &
                           mph_rho1,mph_rho2                , &
                           mph_sten                         , &
                           blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                           blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                           del(DIR_Z)                       , &
                           blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
                           facezData(RH1F_FACE_VAR,:,:,:)   , &
                           facezData(RH2F_FACE_VAR,:,:,:)   , &
                           facezData(SIGM_FACE_VAR,:,:,:) ,&
                           solnData(MDOT_VAR,:,:,:))

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

#ifdef NUCLEATE_BOILING
  do lb = 1,blockCount

     blockID = blockList(lb)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blockID,solnData,CENTER)

     mph_bcFlag(:,:,blockID) = abs(solnData(PFUN_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),blkLimits(LOW,JAXIS),blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  
  end do
#endif

end if

end subroutine
