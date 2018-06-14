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
                       ins_fluxfix_p,&
                   ins_computeQinout,&
                   ins_rescaleVelout,&
                   ins_convectVelout,&
              ins_setInterpValsGcell,&
                           ins_rhs3d_VD,&
                           ins_rhs2d_VD,&
                           ins_rhs2d_PC,&
                           ins_rhs3d_PC,&
                           ins_rhs2d_weno3,&
                           ins_rhs3d_weno3,&
                       ins_predictor_VD,&
                      ins_divergence_PC,&
                       ins_corrector_VD
                     !ins_fluxfix_phi,&


  use IncompNS_data, ONLY : ins_isgs, ins_invRe, ins_intschm, ins_prescoeff, ins_meshMe,&
                            ins_restart, ins_nstep, ins_Qin, ins_Qout, ins_predcorrflg, &
                            ins_convvel, ins_alf, ins_gam, ins_rho, ins_gama, ins_alfa, &
                            ins_rhoa, AB2_SCHM, RK3_SCHM, ins_outflowgridChanged, ins_tlevel, &
                            ins_gravX, ins_gravY,ins_gravZ

  use Grid_Data, ONLY : gr_domainBC 

  use Multiphase_data, only: mph_rho1,mph_rho2,mph_sten,mph_crmx,mph_crmn, &
                             mph_vis1,mph_vis2,mph_lsit,mph_inls,mph_thco1, &
                             mph_thco2,mph_cp1,mph_cp2

   
  use Timers_interface, ONLY : Timers_start, Timers_stop

  use ImBound_interface, ONLY : ImBound

  use Driver_data, ONLY : dr_nstep

  use tree, only : grid_changed, lrefine, surr_blks,nodetype

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use Simulation_data, only: sim_xMin, sim_yMin, sim_zMin
 
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

  real, dimension(GRID_IHI_GC+1,GRID_JHI_GC,GRID_KHI_GC) :: newu
  real, dimension(GRID_IHI_GC,GRID_JHI_GC+1,GRID_KHI_GC) :: newv
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC+1) :: neww

  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_u
  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_v
  real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: flxint_w
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: new_div

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

  integer :: mph_iteration

  real :: xcell,ycell,zcell,xcellX,ycellX,xcellY,ycellY

  character(28) :: filename

  !kpd - for density matching
  integer :: nodetype_perm(MAXBLOCKS)
  integer, save :: mgrid_solveLevelKPD
  integer :: faces(2,MDIM),onBoundary(2,MDIM)

  real :: totaldiv
  real :: pi = 3.1415926535897932384d0
  real :: refPresVal, refPresVal_local

  do lb = 1,blockCount

   blockID = blockList(lb)

   call Grid_getBlkBoundBox(blockId,boundBox)
   bsize(:) = boundBox(2,:) - boundBox(1,:)

   call Grid_getBlkCenterCoords(blockId,coord)

   ! Get blocks dx, dy ,dz:
   call Grid_getDeltas(blockID,del)

   blockID = blockList(lb)
   call Grid_getBlkPtr(blockID,solnData,CENTER)
   call Grid_getBlkPtr(blockID,facexData,FACEX)
   call Grid_getBlkPtr(blockID,faceyData,FACEY)
   call Grid_getBlkPtr(blockID,facezData,FACEZ)

   call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,CENTER)
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

           !zcell = 0.0

           solnData(RTES_VAR,i,j,k) = -12*pi*pi*cos(2*pi*xcell)*cos(2*pi*ycell)*cos(2*pi*zcell)
           solnData(PRES_VAR,i,j,k) = cos(2*pi*xcell)*cos(2*pi*ycell)*cos(2*pi*zcell) - 1.0

        enddo
     enddo
   enddo

   facexData(RH1F_FACE_VAR,:,:,:) = 0.0
   faceyData(RH1F_FACE_VAR,:,:,:) = 0.0
   facezData(RH1F_FACE_VAR,:,:,:) = 0.0
   facexData(RH2F_FACE_VAR,:,:,:) = 1.0
   faceyData(RH2F_FACE_VAR,:,:,:) = 1.0
   facezData(RH2F_FACE_VAR,:,:,:) = 1.0
   solnData(PTES_VAR,:,:,:) = 0.0
   solnData(DELP_VAR,:,:,:) = 0.0

   call Grid_releaseBlkPtr(blockID,solnData,CENTER)
   call Grid_releaseBlkPtr(blockID,solnData,FACEX)
   call Grid_releaseBlkPtr(blockID,solnData,FACEY)
   call Grid_releaseBlkPtr(blockID,solnData,FACEZ)

  end do

  gcMask = .FALSE.
  gcMask(PRES_VAR) = .TRUE.

  gcMask(NUNK_VARS+RH1F_FACE_VAR) = .TRUE.
  gcMask(NUNK_VARS+1*NFACE_VARS+RH1F_FACE_VAR) = .TRUE.
  gcMask(NUNK_VARS+2*NFACE_VARS+RH1F_FACE_VAR) = .TRUE.

  gcMask(NUNK_VARS+RH2F_FACE_VAR) = .TRUE.
  gcMask(NUNK_VARS+1*NFACE_VARS+RH2F_FACE_VAR) = .TRUE.
  gcMask(NUNK_VARS+2*NFACE_VARS+RH2F_FACE_VAR) = .TRUE.

  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,  &
      maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask, &
      selectBlockType=ACTIVE_BLKS)

  newu = 0.
  newv = 0.
  neww = 0.
  flxint_u = 0.
  flxint_v = 0.
  flxint_w = 0.
  new_div = 0.

  nxc = NXB + NGUARD + 1
  nyc = NYB + NGUARD + 1
  nzc = NZB + NGUARD + 1

  do idimn = 1,NDIM
  do ibound = LOW, HIGH

     eachBoundary = 2*(idimn-1)+ibound
     bc_types(eachBoundary) = GRID_PDE_BND_NEUMANN

  enddo
  enddo

  ! SOLUTION OF POISSON EQUATION FOR PRESSURE:
  ! -------- -- ------- -------- --- --------
  call cpu_time(t_startP)
  poisfact = 1.0 
  call Grid_solvePoisson (DELP_VAR, RTES_VAR, bc_types, bc_values, poisfact) 
  call cpu_time(t_stopP)

  if (ins_meshMe .eq. 0) print*,"Total Poisson Solve Time: :",t_stopP-t_startP

  call ins_fluxfix_p(NGUARD,nxc,nyc,nzc,nxc-1,nyc-1,nzc-1,&
                     DELP_VAR,blockCount,blockList)

  gcMask = .FALSE.
  gcMask(DELP_VAR) = .TRUE.
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,  &
      maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask, &
      selectBlockType=ACTIVE_BLKS)

 do lb = 1,blockCount

     blockID = blockList(lb)
     ! Get Index Limits:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     solnData(PTES_VAR,:,:,:) = solnData(PRES_VAR,:,:,:) - solnData(DELP_VAR,:,:,:)
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
 
enddo

  call ins_fluxfix_p(NGUARD,nxc,nyc,nzc,nxc-1,nyc-1,nzc-1,&
                     PTES_VAR,blockCount,blockList)

  gcMask = .FALSE.
  gcMask(PTES_VAR) = .TRUE.
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,  &
      maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask, &
      selectBlockType=ACTIVE_BLKS)

  gcMask(NUNK_VARS+RH1F_FACE_VAR) = .TRUE.
  gcMask(NUNK_VARS+1*NFACE_VARS+RH1F_FACE_VAR) = .TRUE.
  gcMask(NUNK_VARS+2*NFACE_VARS+RH1F_FACE_VAR) = .TRUE.

  gcMask(NUNK_VARS+RH2F_FACE_VAR) = .TRUE.
  gcMask(NUNK_VARS+1*NFACE_VARS+RH2F_FACE_VAR) = .TRUE.
  gcMask(NUNK_VARS+2*NFACE_VARS+RH2F_FACE_VAR) = .TRUE.

END SUBROUTINE ins_ab2rk3_VD
