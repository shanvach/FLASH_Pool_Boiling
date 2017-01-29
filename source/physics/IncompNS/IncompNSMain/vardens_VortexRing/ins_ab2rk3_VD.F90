!!****if* source/physics/IncompNS/IncompNSMain/vardens/ins_ab2rk3_VD
!!
!!
!! NAME
!!
!!  ins_ab2rk3
!!
!!
!! SYNOPSIS
!!
!!  ins_ab2rk3(integer(IN) :: blockCount,
!!             integer(IN) :: blockList(blockCount)
!!             real(IN)    :: timeEndAdv
!!             real(IN)    :: dt)
!!
!!
!! DESCRIPTION
!!
!!  Performs a second order Adams Bashforth or third order Runge-
!!  Kutta step on a fractional step time discretization of the 
!!  Incompressible Navier Stokes flow problem.
!!
!!  The blockList and blockCount arguments tell this routine on
!!  which blocks and on how many to operate.  blockList is an
!!  integer array of size blockCount that contains the local
!!  block numbers of blocks on which to advance.
!!
!!  dt gives the timestep through which this update should advance.
!!
!! ARGUMENTS
!!
!!  blockCount - the number of blocks in blockList
!!  blockList  - array holding local IDs of blocks on which to advance
!!  timeEndAdv - time level at the end of step
!!  dt         - timestep
!!
!!***

subroutine ins_ab2rk3_VD( blockCount, blockList, timeEndAdv, dt)

#include "Flash.h"
#include "ImBound.h"
  ! Modules Use:
#ifdef FLASH_GRID_PARAMESH
  use physicaldata, ONLY : force_consistency,        &
                           interp_mask_unk_res,      &
                           interp_mask_facex_res,    &
                           interp_mask_facey_res,    &
                           interp_mask_facez_res,    &
                           interp_mask_unk,      &
                           interp_mask_facex,    &
                           interp_mask_facey,    &
                           interp_mask_facez
  use workspace, ONLY :    interp_mask_work                           
#endif    

  use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
                             Grid_getListOfBlocks, &
                             Grid_getDeltas,         &
                             Grid_getBlkBC,          &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkIndexLimits, &
                             Grid_fillGuardCells,    &
                             Grid_putFluxData,       &
                             Grid_getFluxData,       &
                             Grid_conserveFluxes,    &
                             Grid_conserveField,     &
                             Grid_updateRefinement,  &
                             Grid_solvePoisson, Grid_getBlkBoundBox, Grid_getBlkCenterCoords

  use gr_interface, ONLY : gr_findMean, gr_findAllNeghID

  use ins_interface, only  :  ins_vt,&
                         ins_vt_WALE,&
                           ins_rhs3d,&
                           ins_rhs2d,&
                       ins_predictor,&
                      ins_divergence,&
                       ins_corrector,&
                         ins_fluxfix,&
                         ins_fluxfixRho1,&
                         ins_fluxfixRho2,&
                         ins_fluxfixSigm,&
                       ins_fluxfix_p,&
                   ins_computeQinout,&
                   ins_rescaleVelout,&
                   ins_convectVelout,&
              ins_setInterpValsGcell,&
                           ins_rhs3d_VD,&
                           ins_rhs2d_VD,&
                       ins_predictor_VD,&
                      ins_divergence_VD,&
                       ins_corrector_VD
                     !ins_fluxfix_phi,&


  use IncompNS_data, ONLY : ins_isgs, ins_invRe, ins_intschm, ins_prescoeff, ins_meshMe,&
                            ins_restart, ins_nstep, ins_Qin, ins_Qout, ins_predcorrflg, &
                            ins_convvel, ins_alf, ins_gam, ins_rho, ins_gama, ins_alfa, &
                            ins_rhoa, AB2_SCHM, RK3_SCHM, ins_outflowgridChanged, ins_tlevel, &
                            ins_gravX, ins_gravY, ins_gravZ

  use Grid_Data, ONLY : gr_domainBC 

  use Multiphase_data, only: mph_rho1,mph_rho2,mph_sten,mph_crmx,mph_crmn, &
                             mph_vis1,mph_vis2,mph_lsit, mph_inls

  use mph_interface, only : mph_KPDcurvature2DAB, mph_KPDcurvature2DC, &
                            mph_KPDadvectWENO3, mph_KPDlsRedistance,  &
                            mph_KPDcurvature3DAB, mph_KPDcurvature3DC,&
                            mph_KPDadvectWENO3_3D, mph_KPDlsRedistance_3D
    
  use Timers_interface, ONLY : Timers_start, Timers_stop

  use ImBound_interface, ONLY : ImBound

  use Driver_data, ONLY : dr_nstep, dr_simTime

  use tree, only : grid_changed, lrefine, surr_blks,nodetype

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

 
  implicit none

#include "constants.h"
#include "IncompNS.h"
!#ifdef FLASH_GRID_PARAMESH
!#include "Multigrid.h"
!#endif
  include "Flash_mpi.h"


  !! ---- Argument List ----------------------------------
  integer, INTENT(INOUT) :: blockCount
  integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList !blockCount
  real,    INTENT(IN) :: timeEndAdv,dt
  !! -----------------------------------------------------

  integer :: temp_grid_changed

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  real, dimension(2,MDIM) :: boundBox


  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)
            
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

  integer :: lb,blockID,ii,jj,kk,ierr,i,j,k,iii,jjj,kkk

  real, dimension(GRID_IHI_GC+1,GRID_JHI_GC,GRID_KHI_GC) :: newu
  real, dimension(GRID_IHI_GC,GRID_JHI_GC+1,GRID_KHI_GC) :: newv
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC+1) :: neww

  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_u
  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_v
  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_w

  integer :: sx,sy,sz,ex,ey,ez

  real dtdxdz,dtdydz,dtdxdy
  
  integer TA(2),count_rate
  real*8  ET

  integer TAIB(2),count_rateIB
  real*8  ETIB

  real maxfp,minfp,maxflb,minflb

  real bsize(MDIM),coord(MDIM)
  integer datasize(MDIM)

  integer nxc, nyc, nzc
  real del(MDIM)

  integer, dimension(6) :: bc_types
  integer :: idimn,ibound,eachBoundary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real, dimension(2,6)  :: bc_values = 0.
  real poisfact,alfadt

  integer ist,itmx
     
  logical, save :: firstcall = .true.

  ! debug VAR:
  integer aa,bb,cc
  real meanPres,meanVelx,meanVely,meanVelz,mndivv,mxdivv,mndivvaux,mxdivvaux

  logical :: gridChanged 
! --------------------------------------------------------------------------
! --------------------------------------------------------------------------
! --------------------------------------------------------------------------
  !kpd
  real :: lsDT,lsT,minCellDiag
  real :: volSum,volSumAll
  real :: rTest1,rTest2

  !- kpd - For Overall Solver Timer... 
  real :: t_startAll,t_stopAll,t_startMP1,t_stopMP1,t_startMP2,t_startMP2a,t_stopMP2, &
          t_startPred,t_stopPred,t_startCorr,t_stopCorr,t_startP,t_stopP

  !- kpd - For Poisson Timer... 
  real :: t_start,t_stop
  integer :: t

  integer :: listofBlocks(MAXBLOCKS)
  integer :: count
  integer :: intval,iOutPress 

  real :: xcell,ycell,zcell,xcellX,ycellX,xcellY,ycellY

  character(28) :: filename

  !kpd - for density matching
  integer :: nodetype_perm(MAXBLOCKS)
  integer, save :: mgrid_solveLevelKPD
  integer :: faces(2,MDIM),onBoundary(2,MDIM)

  real :: Ro, Amp, gamma, xo, yo, zo, Cr, beta, By, ypr, Fy, rii, h1r, Fvort
  real :: Ftime, tpr, alpha, to, tau

! --------------------------------------------------------------------------
! --------------------------------------------------------------------------
! --------------------------------------------------------------------------

  CALL SYSTEM_CLOCK(TA(1),count_rate)  

  newu = 0.
  newv = 0.
  neww = 0.
  flxint_u = 0.
  flxint_v = 0.
  flxint_w = 0.


  nxc = NXB + NGUARD + 1
  nyc = NYB + NGUARD + 1
  nzc = NZB + NGUARD + 1

!!$  write(*,*) 'NXYZB=',NXB,NYB,NZB,NGUARD  
!!$  write(*,*) 'interp_mask_unk=',interp_mask_unk
!!$  write(*,*) 'interp_mask_facex=',interp_mask_facex
!!$  write(*,*) 'interp_mask_facey=',interp_mask_facey
!!$  write(*,*) 'interp_mask_facez=',interp_mask_facez
!!$  write(*,*) 'interp_mask_unk_res=',interp_mask_unk_res
!!$  write(*,*) 'interp_mask_facex_res=',interp_mask_facex_res
!!$  write(*,*) 'interp_mask_facey_res=',interp_mask_facey_res
!!$  write(*,*) 'interp_mask_facez_res=',interp_mask_facez_res


  do idimn = 1,NDIM
  do ibound = LOW, HIGH
     eachBoundary = 2*(idimn-1)+ibound
     select case (gr_domainBC(ibound,idimn))
     case (PERIODIC)
#ifdef FLASH_GRID_UG
        bc_types(eachBoundary) = PERIODIC
#else
        bc_types(eachBoundary) = GRID_PDE_BND_PERIODIC !MG_BND_PERIODIC
#endif
     case (SLIP_INS,NOSLIP_INS,INFLOW_INS,NEUMANN_INS,MOVLID_INS,MOVLIDREV_INS,OUTFLOW_INS)
#ifdef FLASH_GRID_UG
        bc_types(eachBoundary) = OUTFLOW
#else
        bc_types(eachBoundary) = GRID_PDE_BND_NEUMANN !MG_BND_NEUMANN
#endif
     case default
     if (ins_meshMe .eq. MASTER_PE) then
        write(*,*) 'ins_ab2rk3 Error: Boundary Conditions match for Poisson Solver not defined.'
        write(*,*) 'ins_ab2rk3 Error: LOW-HIGH,AXIS=',ibound,idimn
        write(*,*) 'ins_ab2rk3 Error: gr_domainBC(ibound,idimn) =',gr_domainBC(ibound,idimn)
     endif
     call Driver_abortFlash('ins_ab2rk3 Error: BCs do not have matching Poisson solver BCs')
     end select
  enddo
  enddo
 

  ! Select Euler step (for starting) of Adams-Bashforth coefficients
  ! 2nd order Adams Bashforth coefficients:
  if (ins_intschm .eq. AB2_SCHM) then
     ins_gam(1) = 1.5
     ins_gam(2) = 0.0
     ins_gam(3) = 0.0
     ins_rho(1) = -0.5
     ins_rho(2) =  0.0
     ins_rho(3) =  0.0

     itmx = 1
  ! 3rd order Runge Kutta coefficients
  elseif (ins_intschm .eq. RK3_SCHM) then
     ins_gam(1) = 8./15.
     ins_gam(2) = 5./12.
     ins_gam(3) = 3./4.
     ins_rho(1) = 0.0
     ins_rho(2) = -17./60.
     ins_rho(3) = -5./12.

     itmx = 3 
  else
     if (ins_meshMe .eq. MASTER_PE) then
        write(*,*) 'Unknown Incompressible Flow integrator scheme:'
        write(*,*) 'ins_schm=',ins_intschm
        write(*,*) 'where ins_schm=2 Adams-Bashforth, ins_schm=3 Runge-Kutta'
     endif
  endif

  ! Euler coefficients (starting from scratch for Adams-Bashforth):
  if ((ins_nstep .eq. 1) .and. (ins_restart .eqv. .false.) .and. &
      (ins_intschm .eq. AB2_SCHM)) then
     ins_gam(1) = 1.0; ins_gam(2) = 0.0; ins_gam(3) = 0.0
     ins_rho(1) = 0.0; ins_rho(2) = 0.0; ins_rho(3) = 0.0
     itmx = 1
  endif

  ins_alf = ins_gam + ins_rho
 
  ! Set Interpolation values for guardcell-filling:
  call ins_setInterpValsGcell(.true.)

  ins_tlevel = timeEndAdv - dt

!-----------------------------------------------------------------------------------------------
!***********************************************************************************************
!***********************************************************************************************
!***********************************************************************************************
!-----------------------------------------------------------------------------------------------

  !---------------
  ! Timestep Loop:
  !---------------
  do ist = 1,itmx

  !kpd - Start Tital INS_AB2RK3 Timer...
  call cpu_time(t_startAll)

  ins_gama = ins_gam(ist)
  ins_rhoa = ins_rho(ist)
  ins_alfa = ins_alf(ist)

  ins_tlevel = ins_tlevel + ins_alfa*dt

!***********************************************************************************************
!***********************************************************************************************
!***********************************************************************************************

  ! These two subroutine calls ar used in case of outflow BCs, only when NEUMANN_INS and
  ! OUTFLOW_INS are present.
  ! Compute inflow volume ratio: (Not computed on NOT_BOUNDARY, NEUMANN_INS, OUTFLOW_INS)
  call ins_computeQinout( blockCount, blockList, .true., ins_Qin)
  
  ! For OUTFLOW_INS condition compute convective velocity
  call ins_convectVelout( blockCount, blockList, ins_convvel)
  !if(ins_meshMe .eq. MASTER_PE) write(*,*) 'After convect',ins_convvel(HIGH,:)  

!***********************************************************************************************
!***********************************************************************************************
!***********************************************************************************************

  !kpd - Screen Output for Loading...
  call Grid_getListOfBlocks(ALL_BLKS,listofBlocks,count)
  if(ins_meshMe .eq. 0) print*,"The Number of Leaf Blocks (on proc0) is: ",blockCount,"Total: ",count

!***********************************************************************************************
!***********************************************************************************************
!***********************************************************************************************

  ! TURBULENT VISCOSITY COMPUTATION:
  ! --------- --------- -----------
#if NDIM == 3
  if (ins_isgs .NE. 0) then
     do lb = 1,blockCount
        blockID = blockList(lb)

        ! Get blocks dx, dy ,dz:
        call Grid_getDeltas(blockID,del)

        ! Get blocks coord and bsize
        ! Bounding box:
        call Grid_getBlkBoundBox(blockId,boundBox)
        bsize(1:NDIM) = boundBox(2,1:NDIM) - boundBox(1,1:NDIM)

        call Grid_getBlkCenterCoords(blockId,coord)

        ! Point to blocks center and face vars:
        call Grid_getBlkPtr(blockID,solnData,CENTER)
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)
        call Grid_getBlkPtr(blockID,facezData,FACEZ)

        ! calculate turbulent viscosity
        call ins_vt(ins_isgs,NGUARD,nxc,nyc,nzc,              &
                    ins_invRe,del(DIR_X),del(DIR_Y),del(DIR_Z),    &
                    coord,bsize,                                   &
                    facexData,&
                    faceyData,&
                    facezData,&
                    solnData)            

        !call ins_vt_WALE(ins_isgs,NGUARD,nxc,nyc,nzc,              &
        !            ins_invRe,del(DIR_X),del(DIR_Y),del(DIR_Z),    &
        !            coord,bsize,                                   &
        !            facexData,&
        !            faceyData,&
        !            facezData,&
        !            solnData)

        ! Release pointers:
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

     enddo
     ! apply BC and fill guardcells for turbulent viscosity
     gcMask = .FALSE.
     gcMask(TVIS_VAR) = .TRUE.                            ! only turbulent viscosity
     call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)             
  endif
#endif

!!***********************************************************************************************
!!***********************************************************************************************
!    !#############################################################################
!    !-----------------------------------------------------------------------------
!    !-kpd - Fill Guardcells for the distance function before curvature is computed 
!    !       This is done for first iteration only (filled at end of time step) 
!    !-----------------------------------------------------------------------------
!    gcMask = .FALSE.
!    gcMask(DFUN_VAR) = .TRUE.
!#ifdef FLASH_GRID_PARAMESH
!    !intval = 1
!    intval = 2
!    interp_mask_unk = intval;   interp_mask_unk_res = intval;
!    interp_mask_work= intval;
!#endif
!
!    call Grid_fillGuardCells(CENTER,ALLDIR,&
!       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask,selectBlockType=ACTIVE_BLKS)
!      !maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask,selectBlockType=LEAF)
!    !-----------------------------------------------------------------------------
!    !#############################################################################
!!***********************************************************************************************
!!***********************************************************************************************
!!***********************************************************************************************

!kpd - Level Set Initialization...
  if (dr_nstep .eq. 1) then

    !#############################################################################
    !-----------------------------------------------------------------------------
    !-kpd - Fill Guardcells for the distance function before curvature is computed 
    !       This is done for first iteration only (filled at end of time step) 
    !-----------------------------------------------------------------------------
    gcMask = .FALSE.
    gcMask(DFUN_VAR) = .TRUE.
#ifdef FLASH_GRID_PARAMESH
    !intval = 1
    intval = 2
    interp_mask_unk = intval;   interp_mask_unk_res = intval;
    interp_mask_work= intval;
#endif

    call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask,selectBlockType=ACTIVE_BLKS)
      !maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask,selectBlockType=LEAF)
    !-----------------------------------------------------------------------------
    !#############################################################################


   !*********************************************************************************************************
   !- kpd - Level Set Distance Function Initialization (if needed) ******************************************
   !*********************************************************************************************************
   do ii = 1,mph_inls

     !------------------------------
     !- kpd - Level set redistancing 
     !------------------------------

     t = dt
     lsT  = 0.0

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
        if (lb .eq. 1 .AND. ins_meshMe .eq. 0) then
           print*,"Level Set Initialization Iteration # ",ii,minCellDiag,lsDT
        end if

        if (ii.eq.1) solnData(AAJUNK_VAR,:,:,:) = solnData(DFUN_VAR,:,:,:)

        call mph_KPDlsRedistance_3D(solnData(DFUN_VAR,:,:,:), &
                          facexData(VELC_FACE_VAR,:,:,:), &
                          faceyData(VELC_FACE_VAR,:,:,:), &
                          facezData(VELC_FACE_VAR,:,:,:), &
                          del(DIR_X),del(DIR_Y),del(DIR_Z),  &
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS), &
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS), &
                          blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS), &
                          solnData(AAJUNK_VAR,:,:,:), lsDT, blockID )

#elif NDIM == 2


        !--------------------------------------------
        ! Call DFUN re-initialization routine for 2D:
        !--------------------------------------------
        lsDT = MIN(10.0*dt,0.001)
        minCellDiag = SQRT(del(DIR_X)**2.+del(DIR_Y)**2.)
        if (lb .eq. 1 .AND. ins_meshMe .eq. 0) then
           print*,"Level Set Initialization Iteration # ",ii,minCellDiag,lsDT
        end if

        if (ii.eq.1) solnData(AAJUNK_VAR,:,:,:) = solnData(DFUN_VAR,:,:,:)

        call mph_KPDlsRedistance(solnData(DFUN_VAR,:,:,:), &
                          facexData(VELC_FACE_VAR,:,:,:),  &
                          faceyData(VELC_FACE_VAR,:,:,:),  &
                          del(DIR_X),del(DIR_Y),  &
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS), &
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS), &
                          solnData(AAJUNK_VAR,:,:,:), lsDT, minCellDiag)

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
    !-kpd - Fill distance function guard cells after each re-initialization to communicate updates
    gcMask = .FALSE.
    gcMask(DFUN_VAR) = .TRUE.
#ifdef FLASH_GRID_PARAMESH
    !intval = 1
    intval = 2
    interp_mask_unk = intval;   interp_mask_unk_res = intval;
    interp_mask_work= intval;
#endif
    call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)
    !*********************************************************************************************************

      lsT = lsT + lsDT

   end do  ! End do: ii=1,inls
           !==================

   end if  ! End if: dr_nstep = 1
           !=====================

!************************************************************************************************************
!************************************************************************************************************
!************************************************************************************************************
!************************************************************************************************************
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
     !if (ins_nstep .eq. 1)print*,"KPD ins_ab2rk3 bID=",blockID,"blockCount=",blockCount,coord(1),coord(2),coord(3)
     if (ins_nstep .eq. 1)print*,"KPD ins_ab2rk3 bID=",blockID,coord(1),coord(2),coord(3)
     !----------------------------------------------------------

#if NDIM == 2

     !----------------------------------------------------------
     !- kpd - Call 2-D curvature Routine:
     !----------------------------------------------------------
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
                           solnData(VISC_VAR,:,:,:), mph_vis1,mph_vis2,blockID )
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
                           solnData(VISC_VAR,:,:,:), mph_vis1,mph_vis2 )
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
  !!!intval = 2
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
                           facezData(SIGM_FACE_VAR,:,:,:)  )


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
#if NDIM == 3
  gcMask(NUNK_VARS+2*NFACE_VARS+RH1F_FACE_VAR) = .TRUE.    ! rho1z
  gcMask(NUNK_VARS+2*NFACE_VARS+RH2F_FACE_VAR) = .TRUE.    ! rho2z
  gcMask(NUNK_VARS+2*NFACE_VARS+SIGM_FACE_VAR) = .TRUE.    ! Momentum Jump Z
#endif
#ifdef FLASH_GRID_PARAMESH
  intval = 1
  !!!intval = 2
  interp_mask_unk = intval;   interp_mask_unk_res = intval;
  interp_mask_work= intval;
  interp_mask_facex = intval; interp_mask_facex_res = intval;
  interp_mask_facey = intval; interp_mask_facey_res = intval;
  interp_mask_facez = intval; interp_mask_facez_res = intval;
#endif

  !print*,"KPD - Filling Density Guard Cells."
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask,selectBlockType=ACTIVE_BLKS)
  !----------------------------------------------------------------------

!!***********************************************************************************************
!!***********************************************************************************************
!!***********************************************************************************************
!*************************************************************************************************
!*************************************************************************************************
!*************************************************************************************************
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    do lb = 1,count
     blockID = listofBlocks(lb)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_getBlkPtr(blockID,facezData,FACEZ) 
#endif

     facexData(MGW9_FACE_VAR,:,:,:) = 1./(facexData(RH1F_FACE_VAR,:,:,:)+facexData(RH2F_FACE_VAR,:,:,:))
     faceyData(MGW9_FACE_VAR,:,:,:) = 1./(faceyData(RH1F_FACE_VAR,:,:,:)+faceyData(RH2F_FACE_VAR,:,:,:))
#if NDIM == 3
     facezData(MGW9_FACE_VAR,:,:,:) = 1./(facezData(RH1F_FACE_VAR,:,:,:)+facezData(RH2F_FACE_VAR,:,:,:))
#endif

if (blockID .eq. 1) print*,"RHO BEFOR 1",facexData(RH1F_FACE_VAR,20,12,1)+facexData(RH2F_FACE_VAR,20,12,1),"-->",facexData(MGW9_FACE_VAR,20,12,1)
if (blockID .eq. 2) print*,"RHO BEFOR 2",facexData(RH1F_FACE_VAR, 4,12,1)+facexData(RH2F_FACE_VAR, 4,12,1),"-->",facexData(MGW9_FACE_VAR, 4,12,1)
if (blockID .eq. 5) print*,"RHO BEFOR 5",facexData(RH1F_FACE_VAR, 4,4,1)+facexData(RH2F_FACE_VAR, 4,4,1),facexData(RH1F_FACE_VAR, 4,5,1)+facexData(RH2F_FACE_VAR, 4,5,1),"-->",facexData(MGW9_FACE_VAR, 4,4,1),facexData(MGW9_FACE_VAR, 4,5,1)

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

   end do

!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
!!***********************************************************************************************
!!***********************************************************************************************
!!***********************************************************************************************

!***********************************************************************************************
!*************************************************************************************************
!*************************************************************************************************
  ! FIX FLUXES FOR RH1F: (Only for AMR grids)
  ! --- ------ --- -----
#ifdef FLASH_GRID_PARAMESH
  ! Fix fluxes at block boundaries
  call ins_fluxfixRho1(NGUARD,nxc,nyc,nzc,nxc-1,nyc-1,nzc-1,&
                   blockCount,blockList)
#endif
!*************************************************************************************************
!*************************************************************************************************
!*************************************************************************************************

!***********************************************************************************************
!*************************************************************************************************
!*************************************************************************************************
  ! FIX FLUXES FOR RH2F: (Only for AMR grids)
  ! --- ------ --- -----
#ifdef FLASH_GRID_PARAMESH
  ! Fix fluxes at block boundaries
  call ins_fluxfixRho2(NGUARD,nxc,nyc,nzc,nxc-1,nyc-1,nzc-1,&
                   blockCount,blockList)
#endif
!*************************************************************************************************
!*************************************************************************************************
!*************************************************************************************************

!*************************************************************************************************
!*************************************************************************************************
!*************************************************************************************************
  ! FIX FLUXES FOR SIGM: (Only for AMR grids)
  ! --- ------ --- -----
#ifdef FLASH_GRID_PARAMESH
  ! Fix fluxes at block boundaries
  call ins_fluxfixSigm(NGUARD,nxc,nyc,nzc,nxc-1,nyc-1,nzc-1,&
                   blockCount,blockList)
#endif
!*************************************************************************************************
!*************************************************************************************************
!*************************************************************************************************
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    do lb = 1,count
     blockID = listofBlocks(lb)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

if (blockID .eq. 1) print*,"RHO AFTER 1",facexData(RH1F_FACE_VAR,20,12,1)+facexData(RH2F_FACE_VAR,20,12,1),"-->",facexData(MGW9_FACE_VAR,20,12,1)
if (blockID .eq. 2) print*,"RHO AFTER 2",facexData(RH1F_FACE_VAR, 4,12,1)+facexData(RH2F_FACE_VAR, 4,12,1),"-->",facexData(MGW9_FACE_VAR, 4,12,1)
if (blockID .eq. 5) print*,"RHO AFTER 5",facexData(RH1F_FACE_VAR, 4,4,1)+facexData(RH2F_FACE_VAR, 4,4,1),facexData(RH1F_FACE_VAR, 4,5,1)+facexData(RH2F_FACE_VAR, 4,5,1),"-->",facexData(MGW9_FACE_VAR, 4,4,1),facexData(MGW9_FACE_VAR, 4,5,1)

!     do j=4,20
!        do i=4,20
!
!           rTest1 = (facexData(SIGM_FACE_VAR,i+1,j,1)-facexData(SIGM_FACE_VAR,i,j,1))/(del(DIR_X)) + &
!                    (faceyData(SIGM_FACE_VAR,i,j+1,1)-faceyData(SIGM_FACE_VAR,i,j,1))/(del(DIR_Y))
!
!           rTest2 = solnData(SIGP_VAR,i,j,1)
!
!           if ( ABS(rTest1-rTest2) .GT. 0.0001  ) then
!              print*,"DIVERGENCE OF SIGM ERROR!!",blockID,i,j,rTest1,rTest2
!           end if
!
!        end do
!     end do

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)

   end do

   do lb = 1,blockCount
     blockID = blockList(lb)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM ==3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
        do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)

           !rTest1 = (facexData(SIGM_FACE_VAR,i+1,j,1)-facexData(SIGM_FACE_VAR,i,j,1))/(del(DIR_X)) + &
           !         (faceyData(SIGM_FACE_VAR,i,j+1,1)-faceyData(SIGM_FACE_VAR,i,j,1))/(del(DIR_Y))
           !
           !rTest2 = solnData(SIGP_VAR,i,j,1)
           !
           !if ( ABS(rTest1-rTest2) .GT. 0.0001  ) then
           !   print*,"DIVERGENCE OF SIGM ERROR!!",blockID,i,j,rTest1,rTest2
           !end if
           !
           !solnData(SIGP_VAR,i,j,1) = rTest1
#if NDIM == 2
              solnData(SIGP_VAR,i,j,k) = &
                     (facexData(SIGM_FACE_VAR,i+1,j,k)-facexData(SIGM_FACE_VAR,i,j,k))/(del(DIR_X)) + &
                     (faceyData(SIGM_FACE_VAR,i,j+1,k)-faceyData(SIGM_FACE_VAR,i,j,k))/(del(DIR_Y))

#elif NDIM == 3
              solnData(SIGP_VAR,i,j,k) = &
                     (facexData(SIGM_FACE_VAR,i+1,j,k)-facexData(SIGM_FACE_VAR,i,j,k))/(del(DIR_X)) + &
                     (faceyData(SIGM_FACE_VAR,i,j+1,k)-faceyData(SIGM_FACE_VAR,i,j,k))/(del(DIR_Y)) + &
                     (facezData(SIGM_FACE_VAR,i,j,k+1)-facezData(SIGM_FACE_VAR,i,j,k))/(del(DIR_Z))
#endif

           end do
        end do
     end do

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM ==3
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

   end do

!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
!!***********************************************************************************************
!!***********************************************************************************************
!!***********************************************************************************************


!!$  CALL SYSTEM_CLOCK(TA(1),count_rate)  

  ! COMPUTE RIGHT HAND SIDE AND PREDICTOR STEP:
  ! ------- ----- ---- ---- --- --------- ----
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

     call Grid_getBlkBoundBox(blockId,boundBox)
     bsize(:) = boundBox(2,:) - boundBox(1,:)
     call Grid_getBlkCenterCoords(blockId,coord)
     call Grid_getDeltas(blockID,del)

     !=========================================================
     !=========================================================
     !=========================================================
     !KPD - Vortex Ring Source...
     !=========================================================
     Ro = 0.5
     Amp = 350.

     alpha = 500.
     to = 0.05
     tau = 0.04

     beta = 100.0
     yo = -3.*Ro
     By = 0.2*Ro

     gamma = 100.
     xo=0.
     zo=0.
     Cr=0.5

     tpr  =ABS( (dr_simTime-to) ) 

     Ftime = 0.5*(1.+tanh(alpha*(tau-tpr)))

     if (lb .eq. 1) print*,"Vortex Ring Forcing Time:",Ftime

     do k=GRID_KLO_GC,GRID_KHI_GC
        do j=GRID_JLO_GC,GRID_JHI_GC
           do i=GRID_ILO_GC,GRID_IHI_GC

           xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
                   real(i - NGUARD - 1)*del(IAXIS) +   &
                   0.5*del(IAXIS)
           ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                   real(j - NGUARD - 1)*del(JAXIS)  +  &
                   0.5*del(JAXIS)
           zcell  = coord(KAXIS) - bsize(KAXIS)/2.0 +  &
                   real(k - NGUARD - 1)*del(KAXIS)  +  &
                   0.5*del(KAXIS)

           ypr = ABS( (ycell-yo ) )

           Fy = 0.5*(1.+TANH(beta*(By-ypr)))

           rii = sqrt((xcell-xo)**2.+(zcell-zo)**2.)

           h1r = 0.5*(1.+tanh(gamma*(Cr-rii)))

           solnData(AAJUNK_VAR,i,j,k) = Amp*Ftime*Fy*h1r

           !if (ABS(solnData(AAJUNK_VAR,i,j,k)) .gt. 1.e-5) then
           !   print*,"VORTEX:",blockID,i,j,k,solnData(AAJUNK_VAR,i,j,k)
           !end if

           end do
        end do
     end do

     !=========================================================
     !=========================================================
     !=========================================================

     !- kpd - For Predictor Step (newu, newv & neww are RHS)
     call ins_rhs3d_VD (  facexData(VELC_FACE_VAR,:,:,:),            &
                       faceyData(VELC_FACE_VAR,:,:,:),            &
                       facezData(VELC_FACE_VAR,:,:,:),            &
                       solnData(TVIS_VAR,:,:,:),                  &
                       ins_invRe,                                 &
                       blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                       blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                       blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
                       del(DIR_X),del(DIR_Y),del(DIR_Z),newu,newv,neww, &
                       solnData(VISC_VAR,:,:,:),                  &
                       facexData(RH1F_FACE_VAR,:,:,:),            &
                       facexData(RH2F_FACE_VAR,:,:,:),            &
                       faceyData(RH1F_FACE_VAR,:,:,:),            &
                       faceyData(RH2F_FACE_VAR,:,:,:),            &
                       facezData(RH1F_FACE_VAR,:,:,:),            &
                       facezData(RH2F_FACE_VAR,:,:,:),            &
                       ins_gravX, ins_gravY, ins_gravZ, solnData(AAJUNK_VAR,:,:,:) )

     !- kpd - I added this, still a ???
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

     if (ABS(del(DIR_X)-del(DIR_Y)) .gt. 1e-6 .OR. &
         ABS(del(DIR_X)-del(DIR_Z)) .gt. 1e-6 .OR. &
         ABS(del(DIR_Y)-del(DIR_Z)) .gt. 1e-6 ) then
         print*,"Del's:",del(DIR_X),del(DIR_Y),del(DIR_Z)
         call Driver_abortFlash('Cell Spacing Not Equal in X,Y,Z... ins_ab2rk3.f90')
     end if

#elif NDIM ==2

     if (ABS(del(DIR_X)-del(DIR_Y)) .gt. 1e-6 ) then
         print*,"Del's:",del(DIR_X),del(DIR_Y)
        call Driver_abortFlash('Cell Spacing Not Equal in X,Y,Z... ins_ab2rk3.f90')
     end if

     call ins_rhs2d_VD(  facexData(VELC_FACE_VAR,:,:,:),            &
                      faceyData(VELC_FACE_VAR,:,:,:),            &
                      ins_invRe,                                 &
                      blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                      blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                      del(DIR_X),del(DIR_Y),newu,newv, &
                      solnData(VISC_VAR,:,:,:), &
                      facexData(RH1F_FACE_VAR,:,:,:),            &
                      facexData(RH2F_FACE_VAR,:,:,:),            &
                      faceyData(RH1F_FACE_VAR,:,:,:),            &
                      faceyData(RH2F_FACE_VAR,:,:,:),            &
                      ins_gravX, ins_gravY )
     
#endif

     call Grid_getBlkPtr(blockID,facezData,FACEZ)

     call ins_predictor_VD(facexData(VELC_FACE_VAR,:,:,:),&
                        faceyData(VELC_FACE_VAR,:,:,:),&
                        facezData(VELC_FACE_VAR,:,:,:),&
                        newu,newv,neww,                &
                        facexData(RHDS_FACE_VAR,:,:,:),&
                        faceyData(RHDS_FACE_VAR,:,:,:),&
                        facezData(RHDS_FACE_VAR,:,:,:),&
                        solnData(PRES_VAR,:,:,:),      &
                        dt,del(DIR_X),del(DIR_Y),del(DIR_Z),      &
            blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
            blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
            blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
            ins_gama,ins_rhoa,ins_alfa )


     ! save RHS for next step
     facexData(RHDS_FACE_VAR,:,:,:) = newu(:,:,:)
     faceyData(RHDS_FACE_VAR,:,:,:) = newv(:,:,:)


     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)


#if NDIM ==3
     facezData(RHDS_FACE_VAR,:,:,:) = neww(:,:,:)
     !call Grid_releaseBlkPtr(blockID,facezData,FACEZ)    !kpd - I took this out if IF 3d
#endif
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)


  enddo  !end do lb = 1,blockCount

!!$   !CALL SYSTEM_CLOCK(TA(2),count_rate)
!!$   !ET=REAL(TA(2)-TA(1),8)/count_rate
!!$   !write(*,*) 'Predictor time =',ET

  !***********************************************************************************************
  ! APPLY BC AND FILL GUARDCELLS FOR INTERMEDIATE VELOCITIES:
  ! ----- -- --- ---- ---------- --- ------------ ----------
  gcMask = .FALSE.
  gcMask(NUNK_VARS+VELC_FACE_VAR) = .TRUE.                 ! ustar
  gcMask(NUNK_VARS+1*NFACE_VARS+VELC_FACE_VAR) = .TRUE.    ! vstar
#if NDIM == 3
  gcMask(NUNK_VARS+2*NFACE_VARS+VELC_FACE_VAR) = .TRUE.    ! wstar
#endif
  intval = 2   !redundant, but done as a reminder of switch back to quadratic interpolation for u*
  ins_predcorrflg = .true.
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           
  !***********************************************************************************************

!***********************************************************************************************
!*************************************************************************************************
!*************************************************************************************************
  ! FIX FLUXES FOR USTAR: (Only for AMR grids)
  ! --- ------ --- -----
#ifdef FLASH_GRID_PARAMESH
  ! Fix fluxes at block boundaries
  call ins_fluxfix(NGUARD,nxc,nyc,nzc,nxc-1,nyc-1,nzc-1,&
                   blockCount,blockList)
#endif
!*************************************************************************************************
!*************************************************************************************************
!*************************************************************************************************

  ! Compute outflow mass volume ratio: (computed on NEUMANN_INS, OUTFLOW_INS)
  call ins_computeQinout( blockCount, blockList, .false., ins_Qout)
  !if (ins_meshMe .eq. 0) write(*,*) 'Qout before ref=',ins_Qout

  ! Rescale Velocities at outflows for overall conservation: 
  call ins_rescaleVelout(  blockCount, blockList, ins_Qin, ins_Qout)

  CALL SYSTEM_CLOCK(TAIB(1),count_rateIB)
  ! Force Immersed Boundaries:
  call ImBound( blockCount, blockList, ins_alfa*dt,FORCE_FLOW)
  CALL SYSTEM_CLOCK(TAIB(2),count_rateIB)
  ETIB=REAL(TAIB(2)-TAIB(1),8)/count_rateIB
  if (ins_meshMe .eq. MASTER_PE)  write(*,*) 'Total IB Time =',ETIB
 
  ! Compute outflow mass volume ratio: (computed on NEUMANN_INS, OUTFLOW_INS)
  call ins_computeQinout( blockCount, blockList, .false., ins_Qout)
  !if (ins_meshMe .eq. 0) write(*,*) 'Qout after ref=',ins_Qout

!*************************************************************************************************
!*************************************************************************************************
!*************************************************************************************************

  ! DIVERGENCE OF USTAR:
  ! ---------- -- -----
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


#if NDIM ==3
     !call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif
     call Grid_getBlkPtr(blockID,facezData,FACEZ)    !kpd - I moved this
     ! compute divergence of intermediate velocities
     call ins_divergence_VD(facexData(VELC_FACE_VAR,:,:,:),&
                         faceyData(VELC_FACE_VAR,:,:,:),&
                         facezData(VELC_FACE_VAR,:,:,:),&
             blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
             blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
             blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
                       del(DIR_X),del(DIR_Y),del(DIR_Z),&
                       solnData(DUST_VAR,:,:,:) )


     ! Poisson RHS source vector
     solnData(DUST_VAR,                                   &
          blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),     &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),     &
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) =   &
     solnData(DUST_VAR,                                   &
          blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),     &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),     &
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))/(dt*ins_alfa) &
     + &
     solnData(SIGP_VAR,                                   &
          blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),     &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),     &
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM ==3
     !call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif            
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)   !kpd - I moved this

  enddo

!***********************************************************************************************
!***********************************************************************************************
!***********************************************************************************************
!KPD File Output for Testing ONLY
!*********************************************************************************************
!
!    open(unit=123, file='densityWrite.txt',ACTION="write",STATUS="unknown")
!
!    do lb = 1,count
!     blockID = listofBlocks(lb)
!
!     ! Point to blocks center and face vars:
!     call Grid_getBlkPtr(blockID,solnData,CENTER)
!     call Grid_getBlkPtr(blockID,facexData,FACEX)
!     call Grid_getBlkPtr(blockID,faceyData,FACEY)
!
!     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     call Grid_getBlkBoundBox(blockId,boundBox)
!     bsize(:) = boundBox(2,:) - boundBox(1,:)
!     call Grid_getBlkCenterCoords(blockId,coord)
!     call Grid_getDeltas(blockID,del)
!#if NDIM ==3
!     do k=4,20
!#else
!     do k=1,1
!#endif
!     do j=4,20
!        do i=4,20
!
!           xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
!                   real(i - NGUARD - 1)*del(IAXIS) +   &
!                   0.5*del(IAXIS)
!
!           ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
!                   real(j - NGUARD - 1)*del(JAXIS)  +  &
!                   0.5*del(JAXIS)
!#if NDIM ==3
!           zcell  = coord(KAXIS) - bsize(KAXIS)/2.0 +  &
!                   real(k - NGUARD - 1)*del(KAXIS)  +  &
!                   0.5*del(KAXIS)
!#else
!           zcell = 0.0
!#endif
!           xcellX = xcell - 0.5*del(IAXIS)
!           ycellX = ycell
!           xcellY = xcell
!           ycellY = ycell - 0.5*del(JAXIS)
!
!           if (lrefine(lb) .eq. 1) then
!           !if (lrefine(lb) .eq. 2) then
!           !if (lrefine(lb) .eq. 3) then
!           write(123,"(4I5,6F15.6)") lb,i,j,k,xcell,ycell,zcell, &
!                                     facexData(RH1F_FACE_VAR,i,j,1)+facexData(RH2F_FACE_VAR,i,j,1), &
!                                     faceyData(RH1F_FACE_VAR,i,j,1)+faceyData(RH2F_FACE_VAR,i,j,1), &
!                                     solnData(DUST_VAR,i,j,1)
!           end if
!
!        end do
!     end do
!     end do
!     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!     ! Release pointers:
!     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
!     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
!     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
!
!   end do
!
!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************

!  call gr_findMean(DUST_VAR,2,.false.,meanPres)
!  if (ins_meshMe .eq. MASTER_PE) write(*,*) 'Mean Div Ustar A=',meanPres*(dt*ins_alfa)


  ! SOLUTION OF POISSON EQUATION FOR PRESSURE:
  ! -------- -- ------- -------- --- --------
  call cpu_time(t_startP)
  poisfact = 1.0 
  call Grid_solvePoisson (DELP_VAR, DUST_VAR, bc_types, bc_values, poisfact) 
  call cpu_time(t_stopP)
if (ins_meshMe .eq. 0) print*,"Total Poisson Solve Time: :",t_stopP-t_startP

!  call gr_findMean(PRES_VAR,2,.false.,meanPres)
!  if (ins_meshMe .eq. MASTER_PE) write(*,*) 'Mean Pressure=',meanPres
!  call gr_findMean(DELP_VAR,2,.false.,meanPres)
!  if (ins_meshMe .eq. MASTER_PE) write(*,*) 'Mean DeltaP=',meanPres

!***********************************************************************************************
!***********************************************************************************************
!***********************************************************************************************
 
   ! FIX FLUXES FOR dDELP/dxi :
#ifdef FLASH_GRID_UG
  ! Don't Fix Fluxes in block Boundaries
  ! Fill Guardcells for DelP: Used in boundary dDelp/dx fluxes:
  gcMask = .FALSE.
  gcMask(DELP_VAR) = .TRUE.
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,  &
      maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask, &
      selectBlockType=ACTIVE_BLKS)
  ! ---------------------------------------------------------------------
#else
  ! fix fluxes at block boundaries
  ! fix dp gradient fluxes at block boundaries
  call ins_fluxfix_p(NGUARD,nxc,nyc,nzc,nxc-1,nyc-1,nzc-1,&
                     DELP_VAR,blockCount,blockList)
#endif

!***********************************************************************************************
!***********************************************************************************************
!***********************************************************************************************
 
  ! CORRECTOR STEP:
  ! --------- ---
  alfadt = ins_alfa*dt
  do lb = 1,blockCount

     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

#if NDIM == 2
     ! Case 2D take depth to be 1:
     del(DIR_Z) = 1. 
#endif

     dtdydz = alfadt*(del(DIR_Y)*del(DIR_Z))**(-1.)   !(dy*dz)**-1.
     dtdxdz = alfadt*(del(DIR_X)*del(DIR_Z))**(-1.)   !(dx*dz)**-1.
     dtdxdy = alfadt*(del(DIR_X)*del(DIR_Y))**(-1.)   !(dx*dy)**-1.

     ! Get Index Limits:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1   

     ! Positions in face arrays where flux vars have been stored
     sx = NGUARD+1
     sy = NGUARD*K2D+1
     sz = NGUARD*K3D+1
     ex = dataSize(DIR_X)-NGUARD
     ey = dataSize(DIR_Y)-NGUARD*K2D
     ez = dataSize(DIR_Z)-NGUARD*K3D

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

#if NDIM ==3
     !call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif
     call Grid_getBlkPtr(blockID,facezData,FACEZ)  !kpd Moved out of IF

#ifdef FLASH_GRID_UG
     ! UNIFORM Grid:
     ! west face
     facexData(VELC_FACE_VAR,sx,sy:ey,sz:ez) =                  &
          facexData(VELC_FACE_VAR,sx,sy:ey,sz:ez) -             &
          alfadt/del(DIR_X)*(solnData(DELP_VAR,sx,sy:ey,sz:ez)- &
                             solnData(DELP_VAR,sx-1,sy:ey,sz:ez)) * &
          ( facexData(RH1F_FACE_VAR,sx,sy:ey,sz:ez) + facexData(RH2F_FACE_VAR,sx,sy:ey,sz:ez) ) + &
          alfadt*facexData(SIGM_FACE_VAR,sx,sy:ey,sz:ez)

     ! east face
     facexData(VELC_FACE_VAR,ex+1,sy:ey,sz:ez) =                  &
          facexData(VELC_FACE_VAR,ex+1,sy:ey,sz:ez) -             &
          alfadt/del(DIR_X)*(solnData(DELP_VAR,ex+1,sy:ey,sz:ez)- &
                             solnData(DELP_VAR,ex,sy:ey,sz:ez)) * &
          ( facexData(RH1F_FACE_VAR,ex+1,sy:ey,sz:ez) + facexData(RH2F_FACE_VAR,ex+1,sy:ey,sz:ez) ) + &
          alfadt*facexData(SIGM_FACE_VAR,ex+1,sy:ey,sz:ez)

     ! south face
     faceyData(VELC_FACE_VAR,sx:ex,sy,sz:ez) =                  &
          faceyData(VELC_FACE_VAR,sx:ex,sy,sz:ez) -             &
          alfadt/del(DIR_Y)*(solnData(DELP_VAR,sx:ex,sy,sz:ez)- &
                             solnData(DELP_VAR,sx:ex,sy-1,sz:ez)) * &
          ( faceyData(RH1F_FACE_VAR,sx:ex,sy,sz:ez) + faceyData(RH2F_FACE_VAR,sx:ex,sy,sz:ez) ) + &
          alfadt*faceyData(SIGM_FACE_VAR,sx:ex,sy,sz:ez)

     ! north face
     faceyData(VELC_FACE_VAR,sx:ex,ey+1,sz:ez) =                  &
          faceyData(VELC_FACE_VAR,sx:ex,ey+1,sz:ez) -             &
          alfadt/del(DIR_Y)*(solnData(DELP_VAR,sx:ex,ey+1,sz:ez)- &
                             solnData(DELP_VAR,sx:ex,ey,sz:ez)) * &
          ( faceyData(RH1F_FACE_VAR,sx:ex,ey+1,sz:ez) + faceyData(RH2F_FACE_VAR,sx:ex,ey+1,sz:ez) ) + &
          alfadt*faceyData(SIGM_FACE_VAR,sx:ex,ey+1,sz:ez)


#if NDIM == 3
     ! front face
     facezData(VELC_FACE_VAR,sx:ex,sy:ey,sz) =                  &
          facezData(VELC_FACE_VAR,sx:ex,sy:ey,sz) -             &
          alfadt/del(DIR_Z)*(solnData(DELP_VAR,sx:ex,sy:ey,sz)- &
                             solnData(DELP_VAR,sx:ex,sy:ey,sz-1)) * &
          ( facezData(RH1F_FACE_VAR,sx:ex,sy:ey,sz) + facezData(RH2F_FACE_VAR,sx:ex,sy:ey,sz) ) + &
          alfadt*facezData(SIGM_FACE_VAR,sx:ex,sy:ey,sz)

     ! back face
     facezData(VELC_FACE_VAR,sx:ex,sy:ey,ez+1) =                  &
          facezData(VELC_FACE_VAR,sx:ex,sy:ey,ez+1) -             &
          alfadt/del(DIR_Z)*(solnData(DELP_VAR,sx:ex,sy:ey,ez+1)- &
                             solnData(DELP_VAR,sx:ex,sy:ey,ez)) * &
          ( facezData(RH1F_FACE_VAR,sx:ex,sy:ey,ez+1) + facezData(RH2F_FACE_VAR,sx:ex,sy:ey,ez+1) ) + &
          alfadt*facezData(SIGM_FACE_VAR,sx:ex,sy:ey,ez+1)
#endif

#else

     ! AMR GRID:
     ! update block boundary velocities using corrected fluxes
     ! X direction:
     call Grid_getFluxData(blockID, IAXIS, &
                           flxint_u, dataSize)


     ! west face (x goes from sx --> sx)
     facexData(VELC_FACE_VAR,sx,sy:ey,sz:ez) =       &       !u = u - dt/(dy*dz)*flxint_u
          facexData(VELC_FACE_VAR,sx,sy:ey,sz:ez) -  &
          dtdydz*flxint_u(VELC_FLUX,sx,sy:ey,sz:ez)  &
        + alfadt*facexData(SIGM_FACE_VAR,sx,sy:ey,sz:ez)

     ! east face (x goes from ex+1 --> ex+1)
     facexData(VELC_FACE_VAR,ex+1,sy:ey,sz:ez) =       &
          facexData(VELC_FACE_VAR,ex+1,sy:ey,sz:ez) -  &
          dtdydz*flxint_u(VELC_FLUX,ex+1,sy:ey,sz:ez)  &
        + alfadt*facexData(SIGM_FACE_VAR,ex+1,sy:ey,sz:ez)


     ! Y direction:
     ! ------------
     call Grid_getFluxData(blockID, JAXIS, &
                           flxint_v, dataSize)

     ! south face
     faceyData(VELC_FACE_VAR,sx:ex,sy,sz:ez) =       &
          faceyData(VELC_FACE_VAR,sx:ex,sy,sz:ez) -  &
          dtdxdz*flxint_v(VELC_FLUX,sx:ex,sy,sz:ez)  &
        + alfadt*faceyData(SIGM_FACE_VAR,sx:ex,sy,sz:ez)

     ! north face
     faceyData(VELC_FACE_VAR,sx:ex,ey+1,sz:ez) =       &
          faceyData(VELC_FACE_VAR,sx:ex,ey+1,sz:ez) -  &
          dtdxdz*flxint_v(VELC_FLUX,sx:ex,ey+1,sz:ez)  &
        + alfadt*faceyData(SIGM_FACE_VAR,sx:ex,ey+1,sz:ez)

#if NDIM == 3

     ! Z direction:
     call Grid_getFluxData(blockID, KAXIS, &
                           flxint_w, dataSize)

     ! front face
     facezData(VELC_FACE_VAR,sx:ex,sy:ey,sz) =       &
          facezData(VELC_FACE_VAR,sx:ex,sy:ey,sz) -  &
          dtdxdy*flxint_w(VELC_FLUX,sx:ex,sy:ey,sz)  &
        + alfadt*facezData(SIGM_FACE_VAR,sx:ex,sy:ey,sz)

     ! back face
     facezData(VELC_FACE_VAR,sx:ex,sy:ey,ez+1) =       &
          facezData(VELC_FACE_VAR,sx:ex,sy:ey,ez+1) -  &
          dtdxdy*flxint_w(VELC_FLUX,sx:ex,sy:ey,ez+1)  &
        + alfadt*facezData(SIGM_FACE_VAR,sx:ex,sy:ey,ez+1)

#endif

#endif
     ! update divergence-free velocities (not on block boundary)
     call ins_corrector_VD( facexData(VELC_FACE_VAR,:,:,:),&
                         faceyData(VELC_FACE_VAR,:,:,:),&
                         facezData(VELC_FACE_VAR,:,:,:),&
                         facexData(SIGM_FACE_VAR,:,:,:),&
                         faceyData(SIGM_FACE_VAR,:,:,:),&
                         facezData(SIGM_FACE_VAR,:,:,:),&
                         solnData(DELP_VAR,:,:,:),&
                         sx,ex,sy,ey,sz,ez,&
                         dt,del(DIR_X),del(DIR_Y),del(DIR_Z),ins_alfa,  &
                         facexData(RH1F_FACE_VAR,:,:,:),            &
                         facexData(RH2F_FACE_VAR,:,:,:),            &
                         faceyData(RH1F_FACE_VAR,:,:,:),            &
                         faceyData(RH2F_FACE_VAR,:,:,:),            &
                         facezData(RH1F_FACE_VAR,:,:,:),            &
                         facezData(RH2F_FACE_VAR,:,:,:) )

     !- kpd - The final pressure update (For pressure correction ONLY!)
     !        When pressure correction method is not used ins_prescoeff should
     !           be =0.0 and PRES_VAR = DELP_VAR
     solnData(PRES_VAR,:,:,:) = ins_prescoeff*solnData(PRES_VAR,:,:,:) + &
                                solnData(DELP_VAR,:,:,:) 


     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     !call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)  !kpd Moved out of IF
  enddo ! End of corrector loop


  !***********************************************************************************************
  ! FILL GUARDCELLS FOR FINAL VELOCITIES AND PRESSURE:
  ! ---- ---------- --- ----- ---------- --- --------
  ! The pressure fill is used to compute distributed forces on
  ! immersed bodies.
  gcMask = .FALSE.
  gcMask(PRES_VAR) = .TRUE.                                ! pressure
  gcMask(NUNK_VARS+VELC_FACE_VAR) = .TRUE.                 ! u
  gcMask(NUNK_VARS+1*NFACE_VARS+VELC_FACE_VAR) = .TRUE.    ! v
#if NDIM == 3
  gcMask(NUNK_VARS+2*NFACE_VARS+VELC_FACE_VAR) = .TRUE.    ! w
#endif
  ins_predcorrflg = .false.
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)         


!***********************************************************************************************
!***********************************************************************************************
!***********************************************************************************************
!!$  ! FIX FLUXES FOR UFINAL: (Only for AMR grids, needed when using 
!!$  ! boundary restriction of order > 1, and force_consistency_at_srl_interfaces=True)
!!$  ! --- ------ --- -----
!!$#ifdef FLASH_GRID_PARAMESH
!!$  if (force_consistency) then
!!$  ! Fix fluxes at block boundaries
!!$  call ins_fluxfix(NGUARD,nxc,nyc,nzc,nxc-1,nyc-1,nzc-1,&
!!$                   blockCount,blockList)
!!$  endif
!!$#endif
!***********************************************************************************************
!***********************************************************************************************
!***********************************************************************************************

    call cpu_time(t_startMP2)

    !------------------------------------------------------------------
    !- kpd - Advect the multiphase distance function using WENO3 scheme
    !------------------------------------------------------------------

    volSum = 0.0
    volSumAll = 0.0

 !do ii=1,2
  do ii=1,1

    do lb = 1,blockCount
     blockID = blockList(lb)
    !do lb = 1,blockList(blockCount)
    ! blockID = lb

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
  
        !-----------------------------------------------------------------
        !Store Phi at previous time step for RK2
        if (ii.eq.1) solnData(AAJUNK_VAR,:,:,:) = solnData(DFUN_VAR,:,:,:)
        !-----------------------------------------------------------------

        !------------------------------------
        !! Call DFUN advection routine for 3D:
        !------------------------------------
        call mph_KPDadvectWENO3_3D(solnData(DFUN_VAR,:,:,:), &
                          facexData(VELC_FACE_VAR,:,:,:), &
                          faceyData(VELC_FACE_VAR,:,:,:), &
                          facezData(VELC_FACE_VAR,:,:,:), &
                          ins_alfa*dt, &
                          del(DIR_X), &
                          del(DIR_Y), &
                          del(DIR_Z), &
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS), &
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS), &
                          blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),blockID)

        !KPD - Compute the Bubble Volume
        do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
              do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
                 if (solnData(DFUN_VAR,i,j,1) .gt. 0) then
                   volSum = volSum + (del(DIR_X) * del(DIR_Y) * del(DIR_Z))
                 end if
              end do
           end do
        end do

#elif NDIM == 2

        !-----------------------------------------------------------------
        !Store Phi at previous time step for RK2
        if (ii.eq.1) solnData(AAJUNK_VAR,:,:,:) = solnData(DFUN_VAR,:,:,:)
        !-----------------------------------------------------------------

        !------------------------------------
        ! Call DFUN advection routine for 2D:
        !------------------------------------
        call mph_KPDadvectWENO3(solnData(DFUN_VAR,:,:,:), &
                          facexData(VELC_FACE_VAR,:,:,:), &
                          faceyData(VELC_FACE_VAR,:,:,:), &
                          ins_alfa*dt, &
                          del(DIR_X), &
                          del(DIR_Y), &
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),blockID)

        !KPD - Compute the Bubble Volume
        do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
              if (solnData(DFUN_VAR,i,j,1) .gt. 0) then
                volSum = volSum + (del(DIR_X) * del(DIR_Y)) 
              end if
           end do
        end do

#endif

     ! Release pointers:
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM ==3
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
    enddo

    call MPI_Allreduce(volSum, volSumAll, 1, FLASH_REAL,&
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    if (ins_meshMe .eq. 0) print*,"----------------------------------------"
    if (ins_meshMe .eq. 0) print*,"Total Liquid Volume: ",volSumAll
    if (ins_meshMe .eq. 0) print*,"----------------------------------------"

    !********************************************************************************************************
    !-kpd - Fill distance function guard cells before re-initialization to communicate updates
    gcMask = .FALSE.
    gcMask(DFUN_VAR) = .TRUE.
    gcMask(AAJUNK_VAR) = .TRUE.
#ifdef FLASH_GRID_PARAMESH
    !intval = 1
    intval = 2
    interp_mask_unk = intval;   interp_mask_unk_res = intval;
    interp_mask_work = intval;
#endif
    call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)
    !********************************************************************************************************

    if (ii.eq.2) then
       ! Point to blocks center and face vars:
       call Grid_getBlkPtr(blockID,solnData,CENTER)

       solnData(DFUN_VAR,:,:,:) = 0.5* (solnData(AAJUNK_VAR,:,:,:) + solnData(DFUN_VAR,:,:,:))       


       ! Release pointers:
       call Grid_releaseBlkPtr(blockID,solnData,CENTER)        
    end if

  end do   !End ii RK loop


    !********************************************************************************************************
    !-kpd - Fill distance function guard cells before re-initialization to communicate updates
    gcMask = .FALSE.
    gcMask(DFUN_VAR) = .TRUE.
#ifdef FLASH_GRID_PARAMESH
    !intval = 1
    intval = 2
    interp_mask_unk = intval;   interp_mask_unk_res = intval;
    interp_mask_work = intval;
#endif
    call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)
    !********************************************************************************************************

!***********************************************************************************************
!***********************************************************************************************
!***********************************************************************************************

!if (ins_meshMe .eq. 0) print*,"Not filling DFUN GCs, usinf FluxFix"
!
!   ! FIX FLUXES FOR d(DFUN)/dxi :
!#ifdef FLASH_GRID_UG
!  ! ---------------------------------------------------------------------
!#else
!  ! fix fluxes at block boundaries
!  ! fix dp gradient fluxes at block boundaries
!  call ins_fluxfix_phi(NGUARD,nxc,nyc,nzc,nxc-1,nyc-1,nzc-1,&
!                       DFUN_VAR,blockCount,blockList)
!#endif

!***********************************************************************************************
!***********************************************************************************************
!***********************************************************************************************


!************************************************************************************************************
!************************************************************************************************************
!************************************************************************************************************

    call cpu_time(t_startMP2a)

   do ii = 1,mph_lsit

     !------------------------------
     !- kpd - Level set redistancing 
     !------------------------------

     !lsDT = dt
     lsT  = 0.0

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
        lsDT = MIN(50.0*dt,0.05)
        !minCellDiag = SQRT(del(DIR_X)**2.+del(DIR_Y)**2.+del(DIR_Z)**2.)
        minCellDiag = SQRT((SQRT(del(DIR_X)**2.+del(DIR_Y)**2.))**2.+del(DIR_Z)**2.)
        if ( lb .eq. 1 .AND. ins_meshMe .eq. 0) then
           print*,"Level Set Initialization Iteration # ",ii,minCellDiag,lsDT
        end if

!if (blockID .eq. 1) then
!   do jjj=1,22
!      print*,"STRIP",jjj,solnData(DFUN_VAR,12,jjj,8)
!   end do
!end if

        if (ii.eq.1) solnData(AAJUNK_VAR,:,:,:) = solnData(DFUN_VAR,:,:,:)

        call mph_KPDlsRedistance_3D(solnData(DFUN_VAR,:,:,:), &
                          facexData(VELC_FACE_VAR,:,:,:), &
                          faceyData(VELC_FACE_VAR,:,:,:), &
                          facezData(VELC_FACE_VAR,:,:,:), &
                          del(DIR_X),del(DIR_Y),del(DIR_Z),  &
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS), &
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS), &
                          blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS), &
                          solnData(AAJUNK_VAR,:,:,:), lsDT, blockID )

#elif NDIM == 2

        !--------------------------------------------
        ! Call DFUN re-initialization routine for 2D:
        !--------------------------------------------
        lsDT = MIN(10.0*dt,0.005)
        minCellDiag = SQRT(del(DIR_X)**2.+del(DIR_Y)**2.)
        if ( lb .eq. 1 .AND. ins_meshMe .eq. 0) then
           print*,"Level Set Initialization Iteration # ",ii,minCellDiag,lsDT
        end if

        if (ii.eq.1) solnData(AAJUNK_VAR,:,:,:) = solnData(DFUN_VAR,:,:,:)

        call mph_KPDlsRedistance(solnData(DFUN_VAR,:,:,:), &
                          facexData(VELC_FACE_VAR,:,:,:),  &
                          faceyData(VELC_FACE_VAR,:,:,:),  &
                          del(DIR_X),del(DIR_Y),  &
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS), &
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS), &
                          solnData(AAJUNK_VAR,:,:,:), lsDT, minCellDiag)

#endif

        ! Release pointers:
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM ==3
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
     enddo

    !*********************************************************************************************************
    !-kpd - Fill distance function guard cells after each re-initialization to communicate updates
    gcMask = .FALSE.
    gcMask(DFUN_VAR) = .TRUE.
#ifdef FLASH_GRID_PARAMESH
    !intval = 1
    intval = 2
    interp_mask_unk = intval;   interp_mask_unk_res = intval;
#endif
    call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)
    !*********************************************************************************************************

!***********************************************************************************************
!***********************************************************************************************
!***********************************************************************************************
!   ! FIX FLUXES FOR d(DFUN)/dxi :
!#ifdef FLASH_GRID_UG
!  ! ---------------------------------------------------------------------
!#else
!  ! fix fluxes at block boundaries
!  ! fix dp gradient fluxes at block boundaries
!  call ins_fluxfix_phi(NGUARD,nxc,nyc,nzc,nxc-1,nyc-1,nzc-1,&
!                       DFUN_VAR,blockCount,blockList)
!#endif
!***********************************************************************************************
!***********************************************************************************************
!***********************************************************************************************

      lsT = lsT + lsDT

   end do  ! End do: ii=1,lsit

   call cpu_time(t_stopMP2)
   if (ins_meshMe .eq. 0) print*,"Total Multiphase Time: ",t_stopMP2-t_startMP2,t_stopMP2-t_startMP2a

   !print*,"Multiphase 2 Solver Time  ",t_stopMP2-t_startMP2



!***********************************************************************************************
!***********************************************************************************************
!***********************************************************************************************

  ! Compute forces on immersed bodies:
  call ImBound( blockCount, blockList, ins_alfa*dt,COMPUTE_FORCES)

!***********************************************************************************************
!***********************************************************************************************
!***********************************************************************************************

   call cpu_time(t_stopAll)

   if(ins_meshMe .eq. 0) print*,"Total INS Solver Time     ",t_stopAll-t_startAll

!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************

!***********************************************************************************************
!***********************************************************************************************
!***********************************************************************************************
!KPD File Output for Testing ONLY
!*********************************************************************************************

iOutPress = 0

if ((mod(ins_nstep,100) .eq. 0) .AND. (iOutPress .eq. 1)) then

    !open(unit=123, file='pressureWrite.txt',ACTION="write",STATUS="unknown")
    write(filename, '("pressureWrite.txt_",i6.6,"_",i3.3)') ins_nstep,ins_meshMe

    !open(unit=123, file='pressureWrite.txt',ACTION="write",STATUS="unknown")
    open(unit=123, file=filename,ACTION="write",STATUS="unknown")

!    do lb = 1,count
!     blockID = listofBlocks(lb)
     do lb = 1,blockCount
        blockID = blockList(lb)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     call Grid_getBlkBoundBox(blockId,boundBox)
     bsize(:) = boundBox(2,:) - boundBox(1,:)
     call Grid_getBlkCenterCoords(blockId,coord)
     call Grid_getDeltas(blockID,del)
#if NDIM ==3
     do k=4,19
#else
     do k=1,1
#endif
     do j=4,19
        do i=4,19

           xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
                   real(i - NGUARD - 1)*del(IAXIS) +   &
                   0.5*del(IAXIS)

           ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                   real(j - NGUARD - 1)*del(JAXIS)  +  &
                   0.5*del(JAXIS)
#if NDIM ==3
           zcell  = coord(KAXIS) - bsize(KAXIS)/2.0 +  &
                   real(k - NGUARD - 1)*del(KAXIS)  +  &
                   0.5*del(KAXIS)
#else
           zcell = 0.0
#endif
           xcellX = xcell - 0.5*del(IAXIS)
           ycellX = ycell
           xcellY = xcell
           ycellY = ycell - 0.5*del(JAXIS)

           !if (lrefine(lb) .eq. 1) then
           !if (lrefine(lb) .eq. 2) then
           !if (lrefine(lb) .eq. 3) then

           !if (ins_nstep .eq. 1) then
              !write(123,"(4I5,4F15.6)") lb,i,j,k,xcell,ycell,zcell, &
              !                          solnData(PRES_VAR,i,j,1)
              write(123,"(3I5,9F16.11)") lb,i,j,xcell,ycell, &
                                        solnData(PRES_VAR,i,j,1), &
                                        xcellX,ycellX,facexData(VELC_FACE_VAR,i,j,1), &
                                        ycellX,ycellY,faceyData(VELC_FACE_VAR,i,j,1)

           !end if

        end do
     end do
     end do
     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)

   end do

  end if

!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************

  enddo ! End  of time substeps loop


  ! Restore Interpolation values for guardcell-filling:
  call ins_setInterpValsGcell(.false.)

  ! ------------------------------------------------------------------------------------------
  ! Check min max divergence:
  ! ----- --- --- ----------
  mxdivv = -10.**(10.)
  mndivv =  10.**(10.)  
  do lb = 1,blockCount

     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

     !if (blockID .eq. 1) then
     !   print*,"DIV 1",facexData(VELC_FACE_VAR,20,12,1),facexData(VELC_FACE_VAR,19,12,1),faceyData(VELC_FACE_VAR,19,12,1),faceyData(VELC_FACE_VAR,19,13,1),(facexData(VELC_FACE_VAR,20,12,1)-facexData(VELC_FACE_VAR,19,12,1))+(faceyData(VELC_FACE_VAR,19,13,1)-faceyData(VELC_FACE_VAR,19,12,1))
     !end if

#if NDIM == 3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
 
  mxdivv = max( mxdivv,maxval( (facexData(VELC_FACE_VAR,NGUARD+2:nxc,NGUARD+1:nyc-1,NGUARD+1:nzc-1) - &
                    facexData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,NGUARD+1:nzc-1))/del(DIR_X) + &
                   (faceyData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+2:nyc,NGUARD+1:nzc-1) - &
                    faceyData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,NGUARD+1:nzc-1))/del(DIR_Y) + &
                   (facezData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,NGUARD+2:nzc) - &
                    facezData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,NGUARD+1:nzc-1))/del(DIR_Z) ))

  mndivv = min( mndivv,minval( (facexData(VELC_FACE_VAR,NGUARD+2:nxc,NGUARD+1:nyc-1,NGUARD+1:nzc-1) - &
                    facexData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,NGUARD+1:nzc-1))/del(DIR_X) + &
                   (faceyData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+2:nyc,NGUARD+1:nzc-1) - &
                    faceyData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,NGUARD+1:nzc-1))/del(DIR_Y) + &
                   (facezData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,NGUARD+2:nzc) - &
                    facezData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,NGUARD+1:nzc-1))/del(DIR_Z) ))

#elif NDIM == 2

  mxdivv = max( mxdivv,maxval( (facexData(VELC_FACE_VAR,NGUARD+2:nxc,NGUARD+1:nyc-1,1) - &
                    facexData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,1))/del(DIR_X) + &
                   (faceyData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+2:nyc,1) - &
                    faceyData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,1))/del(DIR_Y) ))

  mndivv = min( mndivv,minval( (facexData(VELC_FACE_VAR,NGUARD+2:nxc,NGUARD+1:nyc-1,1) - &
                    facexData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,1))/del(DIR_X) + &
                   (faceyData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+2:nyc,1) - &
                    faceyData(VELC_FACE_VAR,NGUARD+1:nxc-1,NGUARD+1:nyc-1,1))/del(DIR_Y) ))


#endif


     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
  enddo

  mxdivvaux = mxdivv
  call MPI_Allreduce(mxdivvaux, mxdivv, 1, FLASH_REAL,&
                     MPI_MAX, MPI_COMM_WORLD, ierr)

  mndivvaux = mndivv
  call MPI_Allreduce(mndivvaux, mndivv, 1, FLASH_REAL,&
                     MPI_MIN, MPI_COMM_WORLD, ierr)
  

  if (ins_meshMe .eq. MASTER_PE) then
     write(*,*) ' '
     write(*,'(A24,2g14.6)') ' Min , Max  Divergence =',mndivv,mxdivv
  endif

  !----------------------------------------------------------------------------------------------------
  


  CALL SYSTEM_CLOCK(TA(2),count_rate)
  ET=REAL(TA(2)-TA(1),8)/count_rate
  if (ins_meshMe .eq. MASTER_PE)  write(*,*) 'Total AB Step Time =',ET
  

END SUBROUTINE ins_ab2rk3_VD


