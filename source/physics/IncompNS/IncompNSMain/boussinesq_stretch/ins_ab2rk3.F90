!!****if* source/physics/IncompNS/IncompNSMain/constdens/ins_ab2rk3
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

subroutine ins_ab2rk3( blockCount, blockList, timeEndAdv, dt)

#include "Flash.h"
#include "ImBound.h"
#include "SolidMechanics.h"

  use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
                             Grid_getListOfBlocks, &
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

  use gr_interface, ONLY : gr_findMean

  use ins_interface, only  :  ins_vt,&
                           ins_rhs3d,&
                           ins_rhs2d,&
                       ins_predictor,&
                      ins_divergence,&
                       ins_corrector,&
                   !ins_computeQinout,&
                   !ins_rescaleVelout,&
                   !ins_convectVelout,&
              !ins_setInterpValsGcell,&
                      ins_UstarStats!,&
                  !ins_pressgradients

  use IncompNS_data, ONLY : ins_isgs, ins_invsqrtRa_Pr,                                 &
                            ins_intschm, ins_prescoeff, ins_meshMe,                     &
                            ins_restart, ins_nstep, ins_Qin, ins_Qout, ins_predcorrflg, &
                            ins_convvel, ins_alf, ins_gam, ins_rho, ins_gama, ins_alfa, &
                            ins_rhoa, ins_outflowgridChanged, ins_tlevel,               &
                            ins_vardt, rkstep, ins_intschm_type

  use Grid_Data, ONLY : gr_domainBC 
    
  use Timers_interface, ONLY : Timers_start, Timers_stop

  use ImBound_interface, ONLY : ImBound
 
  use SolidMechanics_interface, only : SolidMechanics

  use Driver_interface, only : Driver_getNStep

  implicit none

#include "constants.h"
#include "IncompNS.h"
#include "Flash_mpi.h"

  !! ---- Argument List ----------------------------------
  integer, INTENT(INOUT) :: blockCount
  integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList 
  real,    INTENT(IN) :: timeEndAdv,dt
  !! -----------------------------------------------------

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  real, dimension(2,MDIM) :: boundBox

  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)
            
  real, pointer, dimension(:,:,:,:) :: solnData,scratchData,facexData,faceyData,facezData

  integer :: lb,blockID,ii,jj,kk,ierr,i,j,k

  real, dimension(GRID_IHI_GC+1,GRID_JHI_GC,GRID_KHI_GC) :: newu
  real, dimension(GRID_IHI_GC,GRID_JHI_GC+1,GRID_KHI_GC) :: newv
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC+1) :: neww

  real, dimension(GRID_IHI_GC,3,blockCount) :: iMetrics
  real, dimension(GRID_JHI_GC,3,blockCount) :: jMetrics
  real, dimension(GRID_KHI_GC,3,blockCount) :: kMetrics

  integer TA(2),count_rate
  real*8  ET

  integer TAIB(2),count_rateIB
  real*8  ETIB

  real maxfp,minfp,maxflb,minflb

  real bsize(MDIM),coord(MDIM)
  integer datasize(MDIM)

  integer nxc, nyc, nzc

  integer, dimension(6) :: bc_types
  integer :: idimn,ibound,eachBoundary

  real, dimension(2,6)  :: bc_values = 0.
  real poisfact,alfadt

  integer ist,itmx
     
  logical, save :: firstcall = .true.

  integer :: sm_body_converge

  real :: minu,maxu,minv,maxv,minw,maxw,minp,maxp,mint,maxt,mndivv,mxdivv
  real :: vecminaux(6),vecmaxaux(6),vecmin(6),vecmax(6)

  character(len=6) :: IndNStep
  integer :: NStep

  CALL SYSTEM_CLOCK(TA(1),count_rate)  

  newu = 0.
  newv = 0.
  neww = 0.

  nxc = NXB + NGUARD + 1
  nyc = NYB + NGUARD + 1
  nzc = NZB + NGUARD + 1

  do idimn = 1,NDIM
  do ibound = LOW, HIGH
     eachBoundary = 2*(idimn-1)+ibound
     select case (gr_domainBC(ibound,idimn))
     case (PERIODIC)
        bc_types(eachBoundary) = GRID_PDE_BND_PERIODIC 
     case (SLIP_INS,NOSLIP_INS,INFLOW_INS,NEUMANN_INS,MOVLID_INS,OUTFLOW_INS)
        bc_types(eachBoundary) = GRID_PDE_BND_NEUMANN 
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

#if NDIM==2
  bc_types(5) = bc_types(1)
  bc_types(6) = bc_types(2)
#endif

  ! shift timesteps
  do i = -rkstep,-1
     ins_vardt(i) = ins_vardt(i+1)
  end do
  ins_vardt(0) = dt

  ! Get blk cell metrics by direction from Grid Unit 
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

  ! Select Euler step (for starting) of Adams-Bashforth coefficients
  ! 2nd order Adams Bashforth coefficients (for constant timestep only):
  if (ins_intschm .eq. AB2_SCHM) then
     ins_gam(1) = 1.5
     ins_gam(2) = 0.0
     ins_gam(3) = 0.0
     ins_rho(1) = -0.5
     ins_rho(2) =  0.0
     ins_rho(3) =  0.0

     itmx = 1

  ! 2nd Order Adams-Bashforth coefficents for variable timesteps
  elseif (ins_intschm .eq. AB2_SCHM_V ) then
     ins_gam(1) = 1+ins_vardt(0)/(2.*ins_vardt(-1))
     ins_gam(2) = 0.0
     ins_gam(3) = 0.0
     ins_rho(1) = -ins_vardt(0)/(2.*ins_vardt(-1))
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
      (ins_intschm_type .eq. INS_INTSCHM_MULTISTEP)) then
     ins_gam(1) = 1.0; ins_gam(2) = 0.0; ins_gam(3) = 0.0
     ins_rho(1) = 0.0; ins_rho(2) = 0.0; ins_rho(3) = 0.0
     itmx = 1
  endif

  ins_alf = ins_gam + ins_rho
  ins_tlevel = timeEndAdv - dt

  ! Timestep Loop:
  do ist = 1,itmx

  ins_gama = ins_gam(ist)
  ins_rhoa = ins_rho(ist)
  ins_alfa = ins_alf(ist)

  ins_tlevel = ins_tlevel + ins_alfa*dt

  ! These two subroutine calls ar used in case of outflow BCs, only when NEUMANN_INS and
  ! OUTFLOW_INS are present.
  ! Compute inflow volume ratio: (Not computed on NOT_BOUNDARY, NEUMANN_INS, OUTFLOW_INS)
  !call ins_computeQinout( blockCount, blockList, .true., ins_Qin)
  
  ! For OUTFLOW_INS condition compute convective velocity
  !call ins_convectVelout( blockCount, blockList, ins_convvel)
  !if(ins_meshMe .eq. MASTER_PE) write(*,*) 'After convect',ins_convvel(HIGH,:)  

  ! TURBULENT VISCOSITY COMPUTATION:
  ! --------- --------- -----------
!#if NDIM == 3
!  if (ins_isgs .NE. 0) then
!     do lb = 1,blockCount
!        blockID = blockList(lb)
!
!        ! Get blocks coord and bsize
!        ! Bounding box:
!        call Grid_getBlkBoundBox(blockId,boundBox)
!        bsize(1:NDIM) = boundBox(2,1:NDIM) - boundBox(1,1:NDIM)
!
!        call Grid_getBlkCenterCoords(blockId,coord)
!
!        ! Point to blocks center and face vars:
!        call Grid_getBlkPtr(blockID,solnData,CENTER)
!        call Grid_getBlkPtr(blockID,facexData,FACEX)
!        call Grid_getBlkPtr(blockID,faceyData,FACEY)
!        call Grid_getBlkPtr(blockID,facezData,FACEZ)
!
!        ! calculate turbulent viscosity
!         call ins_vt(ins_isgs,NGUARD,nxc,nyc,nzc,                   &
!                     ins_invsqrtRa_Pr,                              &
!                     del(DIR_X),del(DIR_Y),del(DIR_Z),              &
!                     coord,bsize,                                   &
!                     facexData,&
!                     faceyData,&
!                     facezData,&
!                     solnData)            
!
!        ! Release pointers:
!        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
!        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
!        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
!        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
!
!     enddo
!     ! apply BC and fill guardcells for turbulent viscosity
!     gcMask = .FALSE.
!     gcMask(TVIS_VAR) = .TRUE.                            ! only turbulent viscosity
!     call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
!       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)             
!  endif
!#endif

!!$  CALL SYSTEM_CLOCK(TA(1),count_rate)  


  ! Compute forcing pressure gradients if required:
  !call ins_pressgradients(ins_tlevel,ins_alfa*dt)


  call Timers_start("RightHandSide_Predictor")
  ! COMPUTE RIGHT HAND SIDE AND PREDICTOR STEP:
  ! ------- ----- ---- ---- --- --------- ----
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

     ! compute RHS of momentum equation
     call ins_rhs3d (  facexData(VELC_FACE_VAR,:,:,:),            &
                       faceyData(VELC_FACE_VAR,:,:,:),            & 
                       facezData(VELC_FACE_VAR,:,:,:),            &
                       solnData(TEMP_VAR,:,:,:),                  &
                       solnData(TVIS_VAR,:,:,:),                  &
                       ins_invsqrtRa_Pr,                          &
                       blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                       blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                       blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
                       iMetrics(:,:,lb),                          &
                       jMetrics(:,:,lb),                          &
                       kMetrics(:,:,lb),                          &
                       newu,newv,neww )

#elif NDIM ==2
     ! compute RHS of momentum equation
     call ins_rhs2d(  facexData(VELC_FACE_VAR,:,:,:),             &
                      faceyData(VELC_FACE_VAR,:,:,:),             &
                      solnData(TEMP_VAR,:,:,:),                   &
                      ins_invsqrtRa_Pr,                           &
                      blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS), &
                      blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS), &
                      iMetrics(:,:,lb),jMetrics(:,:,lb),          &
                      newu,newv)
     
#endif
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
     ! compute intermediate velocities
     call ins_predictor(facexData(VELC_FACE_VAR,:,:,:),&
                        faceyData(VELC_FACE_VAR,:,:,:),&
                        facezData(VELC_FACE_VAR,:,:,:),&
                        newu,newv,neww,                &
                        facexData(RHDS_FACE_VAR,:,:,:),&
                        faceyData(RHDS_FACE_VAR,:,:,:),&
                        facezData(RHDS_FACE_VAR,:,:,:),&
                        solnData(PRES_VAR,:,:,:), dt,  &
                        iMetrics(:,LEFT_EDGE,lb), &
                        jMetrics(:,LEFT_EDGE,lb), &
                        kMetrics(:,LEFT_EDGE,lb), &
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
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif


  enddo
  call Timers_stop("RightHandSide_Predictor")

!!$   !CALL SYSTEM_CLOCK(TA(2),count_rate)
!!$   !ET=REAL(TA(2)-TA(1),8)/count_rate
!!$   !write(*,*) 'Predictor time =',ET


  call Timers_start("Gcell_IntermVelocs")
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
  call Timers_stop("Gcell_IntermVelocs")

  ! Copy the original ustar and pressure to velo, preo
  do lb=1,blockCount
     blockID = blockList(lb)
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM ==3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
     facezData(VELO_FACE_VAR,:,:,:) = facezData(VELC_FACE_VAR,:,:,:)
#endif
     ! Velocities
     facexData(VELO_FACE_VAR,:,:,:) = facexData(VELC_FACE_VAR,:,:,:)
     faceyData(VELO_FACE_VAR,:,:,:) = faceyData(VELC_FACE_VAR,:,:,:)
     ! Pressure
     solnData(PREO_VAR,:,:,:)  = solnData(PRES_VAR,:,:,:)

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM ==3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif 
  enddo
  
  ! Begin FSI subiterations:
  sm_body_converge=SM_NOTCONVERGED
  do while(sm_body_converge .eq. SM_NOTCONVERGED)

  ! Copy the original ustar and pressure from velo, preo
  do lb=1,blockCount
     blockID = blockList(lb)
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM ==3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
     facezData(VELC_FACE_VAR,:,:,:) = facezData(VELO_FACE_VAR,:,:,:)
#endif
     ! Velocities
     facexData(VELC_FACE_VAR,:,:,:) = facexData(VELO_FACE_VAR,:,:,:)
     faceyData(VELC_FACE_VAR,:,:,:) = faceyData(VELO_FACE_VAR,:,:,:)
     ! Pressure
     solnData(PRES_VAR,:,:,:)  = solnData(PREO_VAR,:,:,:) 

     call Grid_releaseBlkPtr(blockID,solnData,CENTER) 
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM ==3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif 
  enddo

 
  ! Call Solid Mechanics:
  call Timers_start("SolidMechanics.")
  call SolidMechanics(SM_ADVANCE)  
  call Timers_stop("SolidMechanics.")

!  CALL SYSTEM_CLOCK(TAIB(1),count_rateIB)

  ! Force Immersed Boundaries:
  call Timers_start("Immersed Boundaries.")
  call ImBound( blockCount, blockList, ins_alfa*dt,FORCE_FLOW)
  call Timers_stop("Immersed Boundaries.")

  ! Compute inflow volume ratio: (Not computed on NOT_BOUNDARY, NEUMANN_INS, OUTFLOW_INS)
  !call ins_computeQinout( blockCount, blockList, .true., ins_Qin)  ! Shizhao

  ! Compute DivUstar - delta_mass and print to screen, add to ins_Qin:
  call ins_UstarStats( blockCount, blockList, .true., .true.)

  ! Compute outflow mass volume ratio: (computed on NEUMANN_INS, OUTFLOW_INS)
  !call ins_computeQinout( blockCount, blockList, .false., ins_Qout)

  ! Rescale Velocities at outflows for overall conservation: 
  !call ins_rescaleVelout(  blockCount, blockList, ins_Qin, ins_Qout)

  ! DIVERGENCE OF USTAR:
  ! ---------- -- -----
  do lb = 1,blockCount
     blockID = blockList(lb)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC) 

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)

#if NDIM ==3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif


     ! compute divergence of intermediate velocities
     call ins_divergence(facexData(VELC_FACE_VAR,:,:,:),&
                         faceyData(VELC_FACE_VAR,:,:,:),&
                         facezData(VELC_FACE_VAR,:,:,:),&
             blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
             blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
             blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
                         iMetrics(:,CENTER,lb),         &
                         jMetrics(:,CENTER,lb),         &
                         kMetrics(:,CENTER,lb),         &
                         solnData(DUST_VAR,:,:,:) )


     ! Poisson RHS source vector
     solnData(DUST_VAR,                                   &
          blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),     &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),     &
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) =   &
     solnData(DUST_VAR,                                   &
          blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),     &
          blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),     &
          blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))/(dt*ins_alfa)              


     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM ==3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif            
  enddo

  call Timers_start("Grid_solvePoisson")
  ! SOLUTION OF POISSON EQUATION FOR PRESSURE:
  ! -------- -- ------- -------- --- --------
  poisfact = 1.0 
  call Grid_solvePoisson (DELP_VAR, DUST_VAR, bc_types, bc_values, poisfact) 
  call Timers_stop("Grid_solvePoisson")

  call Timers_start("Gcells_DelP")  
  gcMask = .FALSE.
  gcMask(DELP_VAR) = .TRUE.
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,  &
      maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask, &
      selectBlockType=ACTIVE_BLKS)
  call Timers_stop("Gcells_DelP")
 
  call Timers_start("Corrector")
  ! CORRECTOR STEP:
  ! --------- ---
  alfadt = ins_alfa*dt
  do lb = 1,blockCount

     blockID = blockList(lb)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC) 

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM ==3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif

     ! update divergence-free velocities (not on block boundary)
     !write(*,*) 'ready for corrector on ', ins_meshME
     call ins_corrector( facexData(VELC_FACE_VAR,:,:,:),&
                         faceyData(VELC_FACE_VAR,:,:,:),&
                         facezData(VELC_FACE_VAR,:,:,:),&
                         solnData(DELP_VAR,:,:,:),      & 
             blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
             blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
             blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
                         dt,                            & 
                         iMetrics(:,LEFT_EDGE,blockID), &
                         jMetrics(:,LEFT_EDGE,blockID), &
                         kMetrics(:,LEFT_EDGE,blockID), &
                         ins_alfa)

     ! update pressure
     solnData(PRES_VAR,:,:,:) = ins_prescoeff*solnData(PRES_VAR,:,:,:) + &
                                solnData(DELP_VAR,:,:,:) 

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
  enddo ! End of corrector loop
  call Timers_stop("Corrector")

  call Timers_start("Gcell_FinalVelocsP")
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
  call Timers_stop("Gcell_FinalVelocsP")

  call Timers_start("ImmBoundaries_Forces")
  ! Compute forces on immersed bodies:
  call ImBound( blockCount, blockList, ins_alfa*dt,COMPUTE_FORCES)
  call Timers_stop("ImmBoundaries_Forces")

  ! Call Solid Mechanics:
  call Timers_start("SolidMechanics.")
  call SolidMechanics(SM_CHECKCONVERG,convflag_all=sm_body_converge)  
  call Timers_stop("SolidMechanics.")

  enddo ! End of fsi subiterations

  enddo ! End  of time substeps loop

  ! ------------------------------------------------------------------------------------------
  ! Check min max divergence:
  ! ----- --- --- ----------
  mxdivv = -10.**(10.)
  mndivv =  10.**(10.)  
  maxu   = mxdivv; maxv = maxu; maxw = maxu; maxp = maxu; maxt = maxu
  minu   = mndivv; minv = minu; minw = minu; minp = minu; mint = minu
  do lb = 1,blockCount

     blockID = blockList(lb)

     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,scratchData,SCRATCH_CTR)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif

     call ins_divergence(facexData(VELC_FACE_VAR,:,:,:),&
                         faceyData(VELC_FACE_VAR,:,:,:),&
                         facezData(VELC_FACE_VAR,:,:,:),&
             blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
             blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
             blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
                         iMetrics(:,CENTER,blockID),    &
                         jMetrics(:,CENTER,blockID),    &
                         kMetrics(:,CENTER,blockID),    &
             scratchData(DIVV_SCRATCH_CENTER_VAR,:,:,:) )

  mxdivv = max( mxdivv, maxval(scratchData(DIVV_SCRATCH_CENTER_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)) ) 
  mndivv = min( mndivv, minval(scratchData(DIVV_SCRATCH_CENTER_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)) ) 

  maxu = max(maxu,maxval(facexData(VELC_FACE_VAR,GRID_ILO:GRID_IHI+1,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)))
  minu = min(minu,minval(facexData(VELC_FACE_VAR,GRID_ILO:GRID_IHI+1,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)))

  maxv = max(maxv,maxval(faceyData(VELC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI+1,GRID_KLO:GRID_KHI)))
  minv = min(minv,minval(faceyData(VELC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI+1,GRID_KLO:GRID_KHI)))

  maxp = max(maxp,maxval(solnData(PRES_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)))
  minp = min(minp,minval(solnData(PRES_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)))

  maxt = max(maxt,maxval(solnData(TEMP_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)))
  mint = min(mint,minval(solnData(TEMP_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI)))

#if NDIM == 3

  maxw = max(maxw,maxval(facezData(VELC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI+1)))
  minw = min(minw,minval(facezData(VELC_FACE_VAR,GRID_ILO:GRID_IHI,GRID_JLO:GRID_JHI,GRID_KLO:GRID_KHI+1)))

#elif NDIM == 2

  maxw = 0.0 
  minw = 0.0 

#endif

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,scratchData,SCRATCH_CTR)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

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
  vecmaxaux(6) = maxt
  vecminaux(6) = mint
  
  call MPI_Allreduce(vecmaxaux, vecmax, 6, FLASH_REAL,&
                     MPI_MAX, MPI_COMM_WORLD, ierr)

  call MPI_Allreduce(vecminaux, vecmin, 6, FLASH_REAL,&
                     MPI_MIN, MPI_COMM_WORLD, ierr)
  

  if (ins_meshMe .eq. MASTER_PE) then
     write(*,*) ' '
     write(*,'(A24,2g14.6)') ' Min , Max  U =',vecmin(2),vecmax(2) !minu,maxu
     write(*,'(A24,2g14.6)') ' Min , Max  V =',vecmin(3),vecmax(3) !minv,maxv
#if NDIM == 3
     write(*,'(A24,2g14.6)') ' Min , Max  W =',vecmin(4),vecmax(4) !minw,maxw
#endif
     write(*,'(A24,2g14.6)') ' Min , Max  P =',vecmin(5),vecmax(5) !minp,maxp
     write(*,'(A24,2g14.6)') ' Min , Max  T =',vecmin(6),vecmax(6) !mint,maxt
     write(*,'(A24,2g14.6)') ' Min , Max  Divergence =',vecmin(1),vecmax(1) !mndivv,mxdivv
  endif

  !----------------------------------------------------------------------------------------------------

  CALL SYSTEM_CLOCK(TA(2),count_rate)
  ET=REAL(TA(2)-TA(1),8)/count_rate
  if (ins_meshMe .eq. MASTER_PE)  write(*,*) 'Total AB Step Time =',ET
  

END SUBROUTINE ins_ab2rk3

