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
#include "constants.h"

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
                             GRID_PDE_BND_DIRICHLET, &
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
                             Grid_solvePoisson,      &
                             Grid_getBlkBoundBox,    &
                             Grid_getBlkCenterCoords,&
                             Grid_getCellMetrics

  use gr_interface, ONLY : gr_findMean, gr_findAllNeghID

  use sm_iointerface, only: sm_ioWriteSnapshot,sm_ioWriteParticles, &
                            sm_ioWriteStates, sm_iouttotecplot
 
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
                       ins_fluxfix_p,&
                   ins_computeQinout,&
                   ins_rescaleVelout,&
                   ins_convectVelout,&
              ins_setInterpValsGcell,&
                           ins_rhs3d_VD,&
                           ins_rhs2d_VD,&
                           ins_rhs2d_weno3,&
                           ins_rhs3d_weno3,&
                       ins_predictor_VD,&
                       ins_predictor_VD_IB,&
                       ins_firstCorrector_VD_IB, &
                      ins_divergence_VD,&
                       ins_rhsConstCoef,&
                       ins_correctorConstCoef_IB
                     !ins_fluxfix_phi,&


  use IncompNS_data, ONLY : ins_isgs, ins_invRe, ins_intschm, ins_prescoeff, ins_meshMe,&
                            ins_restart, ins_nstep, ins_Qin, ins_Qout, ins_predcorrflg, &
                            ins_convvel, ins_alf, ins_gam, ins_rho, ins_gama, ins_alfa, &
                            ins_rhoa, AB2_SCHM, RK3_SCHM, ins_outflowgridChanged, ins_tlevel, &
                            ins_gravX, ins_gravY,ins_gravZ

  use Grid_Data, ONLY : gr_domainBC 

  !use Multiphase_data, only: mph_rho1,mph_rho2,mph_sten,mph_crmx,mph_crmn, &
  !                           mph_vis1,mph_vis2,mph_lsit,mph_inls,mph_thco1, &
  !                           mph_thco2,mph_cp1,mph_cp2

  use Multiphase_data, only: mph_rho1, mph_rho2, mph_sten, mph_crmx, mph_crmn, &
                             mph_vis1, mph_vis2, mph_lsit, mph_inls
   
  use Timers_interface, ONLY : Timers_start, Timers_stop

  use ImBound_interface, ONLY : ImBound

  !use Driver_data, ONLY : dr_nstep

  use Driver_data, ONLY : dr_globalMe,dr_simtime,dr_dt,dr_nstep,dr_nstep

  !use tree, only : grid_changed, lrefine, surr_blks,nodetype

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use Simulation_data, only: sim_xMax, sim_xMin, sim_zMax, sim_zMin
  
  use ImBound_data, ONLY: ib_temp_flg, ib_vel_flg, ib_dfun_flg

  use IO_data, ONLY: IO_plotFileIntervalStep

  use SolidMechanics_Data, only : sm_meshMe,sm_meshComm,sm_BodyInfo
 
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

  integer :: lb,blockID,ii,jj,kk,ierr,i,j,k

  real, dimension(GRID_IHI_GC+1,GRID_JHI_GC,GRID_KHI_GC,blockCount) :: tmp_rho1x,tmp_rho2x
  real, dimension(GRID_IHI_GC,GRID_JHI_GC+1,GRID_KHI_GC,blockCount) :: tmp_rho1y,tmp_rho2y
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC+1,blockCount) :: tmp_rho1z,tmp_rho2z

  real, dimension(GRID_IHI_GC+1,GRID_JHI_GC,GRID_KHI_GC) :: newu
  real, dimension(GRID_IHI_GC,GRID_JHI_GC+1,GRID_KHI_GC) :: newv
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC+1) :: neww

  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_u
  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_v
  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_w
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: new_div

  real, dimension(GRID_IHI_GC,3,blockCount) :: dx
  real, dimension(GRID_JHI_GC,3,blockCount) :: dy
  real, dimension(GRID_KHI_GC,3,blockCount) :: dz

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

  real :: rhoz

  integer, dimension(6) :: bc_types
  integer :: idimn,ibound,eachBoundary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real, dimension(2,6)  :: bc_values = 0.
  real poisfact,alfadt

  integer ist,itmx
     
  logical, save :: firstcall = .true.

  ! debug VAR:
  integer aa,bb,cc
  real meanPres,meanVelx,meanVely,meanVelz,mndivv,mxdivv !,mndivvaux,mxdivvaux
  real :: minu,maxu,minv,maxv,minw,maxw,minp,maxp
  real :: vecminaux(5),vecmaxaux(5),vecmin(5),vecmax(5)

  logical :: gridChanged 
! --------------------------------------------------------------------------
! --------------------------------------------------------------------------
! --------------------------------------------------------------------------
  !kpd
  real :: lsDT,lsT,minCellDiag
  real :: volSum,volSumAll

  !- kpd - For Overall Solver Timer... 
  real :: t_startAll,t_stopAll,t_startMP1,t_stopMP1,t_startMP2,t_startMP2a,t_stopMP2, &
          t_startPred,t_stopPred,t_startCorr,t_stopCorr,t_startP,t_stopP

  !- kpd - For Poisson Timer... 
  real :: t_start,t_stop
  integer :: t

  integer :: listofBlocks(MAXBLOCKS)
  integer :: count
  integer :: intval,iOutPress 

  integer :: count1, firstfileflag

  logical :: tecplot_flg

  integer :: mph_iteration

  real :: xcell,ycell,zcell,xcellX,ycellX,xcellY,ycellY

  character(28) :: filename

  !kpd - for density matching
  integer :: nodetype_perm(MAXBLOCKS)
  integer, save :: mgrid_solveLevelKPD
  integer :: faces(2,MDIM),onBoundary(2,MDIM)

  real :: totaldiv
  real :: pi = 3.1415926535897932384d0
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
  new_div = 0.

  tmp_rho1x = 0.
  tmp_rho2x = 0.
  tmp_rho1y = 0.
  tmp_rho2y = 0.
  tmp_rho1z = 0.
  tmp_rho2z = 0.

  nxc = NXB + NGUARD + 1
  nyc = NYB + NGUARD + 1
  nzc = NZB + NGUARD + 1


  ! Get blk cell metrics by direction from Grid Unit 
  do lb = 1,blockCount
     blockID = blockList(lb)
     
     call Grid_getCellMetrics(IAXIS,blockID,LEFT_EDGE, .true.,dx(:,LEFT_EDGE,lb), GRID_IHI_GC) 
     call Grid_getCellMetrics(IAXIS,blockID,CENTER,    .true.,dx(:,CENTER,lb),    GRID_IHI_GC) 
     call Grid_getCellMetrics(IAXIS,blockID,RIGHT_EDGE,.true.,dx(:,RIGHT_EDGE,lb),GRID_IHI_GC) 

     call Grid_getCellMetrics(JAXIS,blockID,LEFT_EDGE, .true.,dy(:,LEFT_EDGE,lb), GRID_JHI_GC) 
     call Grid_getCellMetrics(JAXIS,blockID,CENTER,    .true.,dy(:,CENTER,lb),    GRID_JHI_GC) 
     call Grid_getCellMetrics(JAXIS,blockID,RIGHT_EDGE,.true.,dy(:,RIGHT_EDGE,lb),GRID_JHI_GC) 


     call Grid_getCellMetrics(KAXIS,blockID,LEFT_EDGE, .true.,dz(:,LEFT_EDGE,lb), GRID_KHI_GC) 
     call Grid_getCellMetrics(KAXIS,blockID,CENTER,    .true.,dz(:,CENTER,lb),    GRID_KHI_GC) 
     call Grid_getCellMetrics(KAXIS,blockID,RIGHT_EDGE,.true.,dz(:,RIGHT_EDGE,lb),GRID_KHI_GC) 

  end do

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
     case (SLIP_INS,NOSLIP_INS,INFLOW_INS,NEUMANN_INS,MOVLID_INS,OUTFLOW_INS)
#ifdef FLASH_GRID_UG
        bc_types(eachBoundary) = OUTFLOW
#else
        bc_types(eachBoundary) = GRID_PDE_BND_NEUMANN !MG_BND_NEUMANN !GRID_PDE_BND_DIRICHLET
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

  gcMask = .FALSE.
  gcMask(NUNK_VARS+VELC_FACE_VAR) = .TRUE.                 ! ustar
  gcMask(NUNK_VARS+1*NFACE_VARS+VELC_FACE_VAR) = .TRUE.    ! vstar
#if NDIM == 3
  gcMask(NUNK_VARS+2*NFACE_VARS+VELC_FACE_VAR) = .TRUE.    ! wstar
#endif
  ins_predcorrflg = .false.
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)

  !---------------
  ! Timestep Loop:
  !---------------
  do mph_iteration = 1,1
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
  if(ins_meshMe .eq. MASTER_PE) write(*,*) 'After Qin',ins_Qin  

  ! For OUTFLOW_INS condition compute convective velocity
  call ins_convectVelout( blockCount, blockList, ins_convvel)
  if(ins_meshMe .eq. MASTER_PE) write(*,*) 'After convect',ins_convvel(HIGH,:)  

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

!#if NDIM == 3
  !if (ins_isgs .NE. 0) then
     !do lb = 1,blockCount
        !blockID = blockList(lb)

        ! Get blocks dx, dy ,dz:
        !call Grid_getDeltas(blockID,del)

        ! Get blocks coord and bsize
        ! Bounding box:
        !call Grid_getBlkBoundBox(blockId,boundBox)
        !bsize(1:NDIM) = boundBox(2,1:NDIM) - boundBox(1,1:NDIM)

        !call Grid_getBlkCenterCoords(blockId,coord)

        ! Point to blocks center and face vars:
        !call Grid_getBlkPtr(blockID,solnData,CENTER)
        !call Grid_getBlkPtr(blockID,facexData,FACEX)
        !call Grid_getBlkPtr(blockID,faceyData,FACEY)
        !call Grid_getBlkPtr(blockID,facezData,FACEZ)

        ! calculate turbulent viscosity
        !call ins_vt(ins_isgs,NGUARD,nxc,nyc,nzc,              &
        !            ins_invRe,del(DIR_X),del(DIR_Y),del(DIR_Z),    &
        !            coord,bsize,                                   &
        !            facexData,&
        !            faceyData,&
        !            facezData,&
        !            solnData)            

        !call ins_vt_WALE(ins_isgs,NGUARD,nxc,nyc,nzc,              &
        !            ins_invRe,del(DIR_X),del(DIR_Y),del(DIR_Z),    &
        !            coord,bsize,                                   &
        !            facexData,&
        !            faceyData,&
        !            facezData,&
        !            solnData)

        ! Release pointers:
        !call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        !call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        !call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        !call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

     !enddo
     ! apply BC and fill guardcells for turbulent viscosity
     !gcMask = .FALSE.
     !gcMask(TVIS_VAR) = .TRUE.                            ! only turbulent viscosity
     !call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
     !  maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)             
  !endif
!#endif

  ! COMPUTE RIGHT HAND SIDE AND PREDICTOR STEP:
  ! ------- ----- ---- ---- --- --------- ----
  do lb = 1,blockCount
     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     ! call Grid_getDeltas(blockID,del)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC) 

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,facezData,FACEZ)

#if NDIM == 3
 
     !- avd - For Predictor Step (newu, newv & neww are RHS)
     call ins_rhs3d_weno3(  facexData(VELC_FACE_VAR,:,:,:),            &
                       faceyData(VELC_FACE_VAR,:,:,:),            &
                       facezData(VELC_FACE_VAR,:,:,:),            &
                       solnData(TVIS_VAR,:,:,:),                  &
                       ins_invRe,                                 &
                       blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                       blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                       blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
                       dx(:,:,blockID),dy(:,:,blockID),dz(:,:,blockID),&
                       newu,newv,neww, &
                       solnData(VISC_VAR,:,:,:),                  &
                       facexData(RH1F_FACE_VAR,:,:,:),            &
                       facexData(RH2F_FACE_VAR,:,:,:),            &
                       faceyData(RH1F_FACE_VAR,:,:,:),            &
                       faceyData(RH2F_FACE_VAR,:,:,:),            &
                       facezData(RH1F_FACE_VAR,:,:,:),            &
                       facezData(RH2F_FACE_VAR,:,:,:),            &
                       ins_gravX, ins_gravY, ins_gravZ)

!     if (ABS(del(DIR_X)-del(DIR_Y)) .gt. 1e-6 .OR. &
!         ABS(del(DIR_X)-del(DIR_Z)) .gt. 1e-6 .OR. &
!         ABS(del(DIR_Y)-del(DIR_Z)) .gt. 1e-6 ) then
!         print*,"Del's:",del(DIR_X),del(DIR_Y),del(DIR_Z)
!         call Driver_abortFlash('Cell Spacing Not Equal in X,Y,Z... ins_ab2rk3.f90')
!     end if

#elif NDIM ==2

!     if (ABS(del(DIR_X)-del(DIR_Y)) .gt. 1e-6 ) then
!         print*,"Del's:",del(DIR_X),del(DIR_Y)
!        call Driver_abortFlash('Cell Spacing Not Equal in X,Y,Z... ins_ab2rk3.f90')
!     end if

     ! ML - GFM for velocity jump condition
     call ins_rhs2d_weno3(  facexData(VELC_FACE_VAR,:,:,:),            &
                      faceyData(VELC_FACE_VAR,:,:,:),            &
                      ins_invRe,                                 &
                      blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                      blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                      dx(:,:,blockID),dy(:,:,blockID),newu,newv, &
                      solnData(VISC_VAR,:,:,:), &
                      facexData(RH1F_FACE_VAR,:,:,:),            &
                      facexData(RH2F_FACE_VAR,:,:,:),            &
                      faceyData(RH1F_FACE_VAR,:,:,:),            &
                      faceyData(RH2F_FACE_VAR,:,:,:),            &
                      ins_gravX, ins_gravY)
     
#endif

     call ins_predictor_VD_IB(facexData(VELC_FACE_VAR,:,:,:),&
                        faceyData(VELC_FACE_VAR,:,:,:),&
                        facezData(VELC_FACE_VAR,:,:,:),&
                        newu,newv,neww,                &
                        facexData(RHDS_FACE_VAR,:,:,:),&
                        faceyData(RHDS_FACE_VAR,:,:,:),&
                        facezData(RHDS_FACE_VAR,:,:,:),&
                        solnData(PRES_VAR,:,:,:),      &
                        dt, &
            blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
            blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
            blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
            ins_gama,ins_rhoa,ins_alfa, &
                        facexData(POLD_FACE_VAR,:,:,:),&
                        faceyData(POLD_FACE_VAR,:,:,:),&
                        facezData(POLD_FACE_VAR,:,:,:))

     ! save RHS for next step
     facexData(RHDS_FACE_VAR,:,:,:) = newu(:,:,:)
     faceyData(RHDS_FACE_VAR,:,:,:) = newv(:,:,:)
#if NDIM ==3
     facezData(RHDS_FACE_VAR,:,:,:) = neww(:,:,:)
#endif

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
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
  ins_predcorrflg = .true.
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           
  !***********************************************************************************************

!***********************************************************************************************
!*************************************************************************************************
!*************************************************************************************************

  ! Compute outflow mass volume ratio: (computed on NEUMANN_INS, OUTFLOW_INS)
  call ins_computeQinout( blockCount, blockList, .false., ins_Qout)
  if (ins_meshMe .eq. 0) write(*,*) 'Qout before ref=',ins_Qout
  !if (ins_meshMe .eq. 0) stop
  ! Rescale Velocities at outflows for overall conservation: 
  call ins_rescaleVelout(  blockCount, blockList, ins_Qin, ins_Qout) !- ML: commented out due to error with neumann_ins bc?

  !ib_temp_flg = .false.
  ib_vel_flg  = .true.
  ib_dfun_flg = .false.

  CALL SYSTEM_CLOCK(TAIB(1),count_rateIB)
  ! Force Immersed Boundaries:
  call ImBound( blockCount, blockList, ins_alfa*dt,FORCE_FLOW) 
  CALL SYSTEM_CLOCK(TAIB(2),count_rateIB)
  ETIB=REAL(TAIB(2)-TAIB(1),8)/count_rateIB
  if (ins_meshMe .eq. MASTER_PE)  write(*,*) 'Total IB Time =',ETIB

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
     call Grid_getBlkPtr(blockID,facezData,FACEZ)

    !if(dr_nstep==1 .and. ins_meshMe == 2) then
    !open(1, file = 'divu2.txt', action = 'write', access = 'append') !status ='replace'
    !    do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
    !       do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
    !         Mdens = (rho1x(i,j,k) + rho2x(i,j,k))
    !         write(1,*) i,j,facexData(VELC_FACE_VAR,i,j,1),ins_meshMe
    !          write(1,*) i,j,facexData(VELC_FACE_VAR,i,j,1),faceyData(VELC_FACE_VAR,i,j,1)
    !                               ((2.0*po(i,j,k)-1.0*poo(i,j,k))-&
    !                                        (2.0*po(i-1,j,k)-1.0*poo(i-1,j,k)))/dx
    !       enddo
    !     enddo
    !    write(1,*) "---"
    ! close(1)              
    ! endif

     call ins_firstCorrector_VD_IB( facexData(VELC_FACE_VAR,:,:,:),&
                         faceyData(VELC_FACE_VAR,:,:,:),&
                         facezData(VELC_FACE_VAR,:,:,:),&
                         blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                         blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                         blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
                         dt, &
                         facexData(POLD_FACE_VAR,:,:,:),            &
                         faceyData(POLD_FACE_VAR,:,:,:),            &
                         facezData(POLD_FACE_VAR,:,:,:))

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
  end do

  gcMask = .FALSE.
  gcMask(NUNK_VARS+VELC_FACE_VAR) = .TRUE.                 ! ustar
  gcMask(NUNK_VARS+1*NFACE_VARS+VELC_FACE_VAR) = .TRUE.    ! vstar
#if NDIM == 3
  gcMask(NUNK_VARS+2*NFACE_VARS+VELC_FACE_VAR) = .TRUE.    ! wstar
#endif
  ins_predcorrflg = .false.
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           

  ! Compute outflow mass volume ratio: (computed on NEUMANN_INS, OUTFLOW_INS)
  call ins_computeQinout( blockCount, blockList, .false., ins_Qout)
  if (ins_meshMe .eq. 0) write(*,*) 'Qout after ref=',ins_Qout
  !if(ins_meshMe .eq. MASTER_PE) stop

!*************************************************************************************************
!*************************************************************************************************
!*************************************************************************************************

  rhoz = min(mph_rho1,mph_rho2)/mph_rho2

  ! DIVERGENCE OF USTAR:
  ! ---------- -- -----
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
     call Grid_getBlkPtr(blockID,facezData,FACEZ)

     ! compute divergence of intermediate velocities

     !-avd-compute divergence for phase problems
     call ins_divergence_VD(facexData(VELC_FACE_VAR,:,:,:),&
                         faceyData(VELC_FACE_VAR,:,:,:),&
                         facezData(VELC_FACE_VAR,:,:,:),&
             blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
             blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
             blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
                       dx(:,CENTER,blockID),dy(:,CENTER,blockID),&
                       dz(:,CENTER,blockID),&
                       solnData(DUST_VAR,:,:,:))
     ! Poisson RHS source vector


     call ins_rhsConstCoef (blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                            blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                            blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
                            dx(:,:,blockID),dy(:,:,blockID),&
                            dz(:,:,blockID),dt,mph_rho1,mph_rho2,&
                            solnData(PRES_VAR,:,:,:),solnData(POOD_VAR,:,:,:),&
                            solnData(SPCD_VAR,:,:,:),&
                            facexData(SIGO_FACE_VAR,:,:,:),&
                            faceyData(SIGO_FACE_VAR,:,:,:),&
                            facezData(SIGO_FACE_VAR,:,:,:),&
                            facexData(SIGOO_FACE_VAR,:,:,:),&
                            faceyData(SIGOO_FACE_VAR,:,:,:),&
                            facezData(SIGOO_FACE_VAR,:,:,:),&
                            facexData(RH1F_FACE_VAR,:,:,:),&
                            facexData(RH2F_FACE_VAR,:,:,:),&
                            faceyData(RH1F_FACE_VAR,:,:,:),&
                            faceyData(RH2F_FACE_VAR,:,:,:),&
                            facezData(RH1F_FACE_VAR,:,:,:),&
                            facezData(RH2F_FACE_VAR,:,:,:),rhoz,&
                            solnData(PFUN_VAR,:,:,:),&
                            solnData(DFUN_VAR,:,:,:),&
                            solnData(DUST_VAR,:,:,:),&
                            tmp_rho1x(:,:,:,blockID),tmp_rho2x(:,:,:,blockID),&
                            tmp_rho1y(:,:,:,blockID),tmp_rho2y(:,:,:,blockID),&
                            tmp_rho1z(:,:,:,blockID),tmp_rho2z(:,:,:,blockID),blockID)

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)            
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  enddo

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
     !call Grid_getDeltas(blockID,del)

!#if NDIM == 2
     ! Case 2D take depth to be 1:
     !del(DIR_Z) = 1. 
!#endif

     ! Get Index Limits:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

      
     datasize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1   

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
 
     call ins_correctorConstCoef_IB(facexData(VELC_FACE_VAR,:,:,:),&
                           faceyData(VELC_FACE_VAR,:,:,:),&
                           facezData(VELC_FACE_VAR,:,:,:),&
                           facexData(SIGC_FACE_VAR,:,:,:),&
                           faceyData(SIGC_FACE_VAR,:,:,:),&
                           facezData(SIGC_FACE_VAR,:,:,:),&
                           facexData(SIGO_FACE_VAR,:,:,:),&
                           faceyData(SIGO_FACE_VAR,:,:,:),&
                           facezData(SIGO_FACE_VAR,:,:,:),&
                           facexData(SIGOO_FACE_VAR,:,:,:),&
                           faceyData(SIGOO_FACE_VAR,:,:,:),&
                           facezData(SIGOO_FACE_VAR,:,:,:),&
                           facexData(RH1F_FACE_VAR,:,:,:),&
                           facexData(RH2F_FACE_VAR,:,:,:),&
                           faceyData(RH1F_FACE_VAR,:,:,:),&
                           faceyData(RH2F_FACE_VAR,:,:,:),&
                           facezData(RH1F_FACE_VAR,:,:,:),&
                           facezData(RH2F_FACE_VAR,:,:,:),&
                           solnData(DELP_VAR,:,:,:),&
                           solnData(PRES_VAR,:,:,:),&
                           solnData(POOD_VAR,:,:,:),&
                           blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                           blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                           blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
                           rhoz,dt,dx(:,LEFT_EDGE,blockID),dy(:,LEFT_EDGE,blockID),&
                           dz(:,LEFT_EDGE,blockID),mph_rho1,mph_rho2,&
                           solnData(PFUN_VAR,:,:,:),&
                           solnData(DFUN_VAR,:,:,:),ins_alfa,&
                           facexData(POLD_FACE_VAR,:,:,:),&
                           faceyData(POLD_FACE_VAR,:,:,:),&
                           facezData(POLD_FACE_VAR,:,:,:),&
                           tmp_rho1x(:,:,:,blockID),tmp_rho2x(:,:,:,blockID),&
                           tmp_rho1y(:,:,:,blockID),tmp_rho2y(:,:,:,blockID),&
                           tmp_rho1z(:,:,:,blockID),tmp_rho2z(:,:,:,blockID),blockID)


     !print*,maxval(facexData(POLD_FACE_VAR,:,:,:))
     !stop
     facexData(SIGOO_FACE_VAR,:,:,:) = facexData(SIGO_FACE_VAR,:,:,:)
     facexData(SIGO_FACE_VAR,:,:,:) = facexData(SIGC_FACE_VAR,:,:,:)
     faceyData(SIGOO_FACE_VAR,:,:,:) = faceyData(SIGO_FACE_VAR,:,:,:)
     faceyData(SIGO_FACE_VAR,:,:,:) = faceyData(SIGC_FACE_VAR,:,:,:)
#if NDIM == 3
     facezData(SIGOO_FACE_VAR,:,:,:) = facezData(SIGO_FACE_VAR,:,:,:)
     facezData(SIGO_FACE_VAR,:,:,:) = facezData(SIGC_FACE_VAR,:,:,:)
#endif
     solnData(POOD_VAR,:,:,:) = solnData(PRES_VAR,:,:,:)  
 
     solnData(PRES_VAR,:,:,:) = ins_prescoeff*solnData(PRES_VAR,:,:,:) + &
                                solnData(DELP_VAR,:,:,:) 
 
     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  enddo ! End of corrector loop


  !***********************************************************************************************
  ! FILL GUARDCELLS FOR FINAL VELOCITIES AND PRESSURE:
  ! ---- ---------- --- ----- ---------- --- --------
  ! The pressure fill is used to compute distributed forces on
  ! immersed bodies.
  gcMask = .FALSE.
  gcMask(PRES_VAR) = .TRUE.                                ! pressure
  gcMask(POOD_VAR) = .TRUE.                                ! pressure
  gcMask(NUNK_VARS+VELC_FACE_VAR) = .TRUE.                 ! u
  gcMask(NUNK_VARS+1*NFACE_VARS+VELC_FACE_VAR) = .TRUE.    ! v
#if NDIM == 3
  gcMask(NUNK_VARS+2*NFACE_VARS+VELC_FACE_VAR) = .TRUE.    ! w
#endif
  gcMask(NUNK_VARS+POLD_FACE_VAR) = .TRUE.                 ! u
  gcMask(NUNK_VARS+1*NFACE_VARS+POLD_FACE_VAR) = .TRUE.    ! v
#if NDIM == 3
  gcMask(NUNK_VARS+2*NFACE_VARS+POLD_FACE_VAR) = .TRUE.    ! w
#endif 
  ins_predcorrflg = .false.
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)         

  call cpu_time(t_stopAll)

   if(ins_meshMe .eq. 0) print*,"Total INS Solver Time     ",t_stopAll-t_startAll

!*********************************************************************************************
!*********************************************************************************************

  enddo ! End  of time substeps loop
  enddo

  CALL SYSTEM_CLOCK(TAIB(1),count_rateIB)
  ! Force Immersed Boundaries:
  call ImBound( blockCount, blockList, ins_alfa*dt,COMPUTE_FORCES) 
  CALL SYSTEM_CLOCK(TAIB(2),count_rateIB)
  ETIB=REAL(TAIB(2)-TAIB(1),8)/count_rateIB
  if (ins_meshMe .eq. MASTER_PE)  write(*,*) 'Total IB Time =',ETIB

  ! Restore Interpolation values for guardcell-filling:
  call ins_setInterpValsGcell(.false.)


 ! ------------------------------------------------------------------------------------------
  ! Check min max divergence:
  ! ----- --- --- ----------
  mxdivv = -10.**(10.)
  mndivv =  10.**(10.)
  maxu   = mxdivv; maxv = maxu; maxw = maxu; maxp = maxu;
  minu   = mndivv; minv = minu; minw = minu; minp = minu;
  do lb = 1,blockCount

     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     !call Grid_getDeltas(blockID,del)

     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,facezData,FACEZ)

     call ins_divergence_VD(facexData(VELC_FACE_VAR,:,:,:),&
                         faceyData(VELC_FACE_VAR,:,:,:),&
                         facezData(VELC_FACE_VAR,:,:,:),&
             blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
             blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
             blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
                       dx(:,CENTER,blockID),dy(:,CENTER,blockID),&
                       dz(:,CENTER,blockID),new_div)

  mxdivv = max( mxdivv,maxval(new_div))

  mndivv = min( mndivv,minval(new_div))

  maxu = max(maxu,maxval(facexData(VELC_FACE_VAR,GRID_ILO:GRID_IHI+1,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)))
  minu = min(minu,minval(facexData(VELC_FACE_VAR,GRID_ILO:GRID_IHI+1,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)))

  maxv = max(maxv,maxval(faceyData(VELC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI+1,GRID_KLO:GRID_KHI)))
  minv = min(minv,minval(faceyData(VELC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI+1,GRID_KLO:GRID_KHI)))

#if NDIM == 3
  maxw = max(maxw,maxval(facezData(VELC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI+1)))
  minw = min(minw,minval(facezData(VELC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI+1)))
#elif NDIM == 2
  maxw = 0.0
  minw = 0.0
#endif

  maxp = max(maxp,maxval(solnData(PRES_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)))
  minp = min(minp,minval(solnData(PRES_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)))

  call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  call Grid_releaseBlkPtr(blockID,facexData,FACEX)
  call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
  call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  enddo

  vecmaxaux(1) = mxdivv
  vecminaux(1) = mndivv
  vecmaxaux(2) = maxu
  vecminaux(2) = minu
  vecmaxaux(3) = maxv
  vecminaux(3) = minv
  vecmaxaux(4) = maxw
  vecminaux(4) = minw
  vecmaxaux(5) = maxp
  vecminaux(5) = minp

  call MPI_Allreduce(vecmaxaux, vecmax, 5, FLASH_REAL,&
                     MPI_MAX, MPI_COMM_WORLD, ierr)
 
  call MPI_Allreduce(vecminaux, vecmin, 5, FLASH_REAL,&
                     MPI_MIN, MPI_COMM_WORLD, ierr)

  if (ins_meshMe .eq. MASTER_PE) then
     write(*,*) ' '
     write(*,'(A24,2g14.6)') ' Min , Max  U =',vecmin(2),vecmax(2) !minu,maxu
     write(*,'(A24,2g14.6)') ' Min , Max  V =',vecmin(3),vecmax(3) !minv,maxv
     write(*,'(A24,2g14.6)') ' Min , Max  W =',vecmin(4),vecmax(4) !minw,maxw
     write(*,'(A24,2g14.6)') ' Min , Max  P =',vecmin(5),vecmax(5) !minp,maxp
     write(*,'(A24,2g14.6)') ' Min , Max  Divergence =',vecmin(1),vecmax(1) !mndivv,mxdivv
  endif

  !----------------------------------------------------------------------------------------------------


  CALL SYSTEM_CLOCK(TA(2),count_rate)
  ET=REAL(TA(2)-TA(1),8)/count_rate
  if (ins_meshMe .eq. MASTER_PE)  write(*,*) 'Total AB Step Time =',ET
  

END SUBROUTINE ins_ab2rk3_VD


