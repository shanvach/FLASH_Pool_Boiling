subroutine ins_MFcorrection( blockCount, blockList,timeEndAdv, dt)

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


  use ins_interface, ONLY: ins_divergence_PC

  use IncompNS_data, ONLY : ins_isgs, ins_invRe, ins_intschm, ins_prescoeff,ins_meshMe,&
                            ins_restart, ins_nstep, ins_Qin, ins_Qout, ins_predcorrflg, &
                            ins_convvel, ins_alf, ins_gam, ins_rho, ins_gama, ins_alfa, &
                            ins_rhoa, AB2_SCHM, RK3_SCHM, ins_outflowgridChanged, ins_tlevel, &
                            ins_gravX, ins_gravY, ins_gravZ

  use Multiphase_data, only: mph_rho1,mph_rho2,mph_sten,mph_crmx,mph_crmn, &
                             mph_vis1,mph_vis2,mph_lsit,mph_inls,mph_thco1, &
                             mph_thco2,mph_cp1,mph_cp2


  use Driver_data, ONLY : dr_nstep

  use Grid_Data, ONLY : gr_domainBC

  implicit none

#include "constants.h"
#include "IncompNS.h"

  include "Flash_mpi.h"

  !! ---- Argument List ----------------------------------
  integer, INTENT(INOUT) ::  blockCount
  integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList
  real,    INTENT(IN) :: dt, timeEndAdv
  !! -----------------------------------------------------

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  real, dimension(2,MDIM) :: boundBox

  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

  integer :: lb,blockID,ii,jj,kk,ierr,i,j,k

  integer :: sx,sy,sz,ex,ey,ez

  integer :: ix1,ix2,jy1,jy2,kz1,kz2

  real dtdxdz,dtdydz,dtdxdy

  real bsize(MDIM),coord(MDIM)

  integer datasize(MDIM)

  integer nxc, nyc, nzc

  real del(MDIM)

  real, dimension(2,6)  :: bc_values = 0.

  real poisfact,alfadt

  integer :: listofBlocks(MAXBLOCKS)
  integer :: count
  integer :: intval,iOutPress

  !kpd - for density matching
  integer :: nodetype_perm(MAXBLOCKS)
  integer, save :: mgrid_solveLevelKPD
  integer :: faces(2,MDIM),onBoundary(2,MDIM)

  integer, dimension(6) :: bc_types
  integer :: idimn,ibound,eachBoundary


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
        bc_types(eachBoundary) = GRID_PDE_BND_DIRICHLET !MG_BND_NEUMANN
!GRID_PDE_BND_DIRICHLET
#endif
     case default
     end select
  enddo
  enddo

   do lb = 1,blockCount

      blockID = blockList(lb)


        ! Get blocks dx, dy ,dz:
        call Grid_getDeltas(blockID,del)

        ! Get blocks coord and bsize
        ! Bounding box:
        call Grid_getBlkBoundBox(blockId,boundBox)
        bsize(1:NDIM) = boundBox(2,1:NDIM) - boundBox(1,1:NDIM)

        call Grid_getBlkCenterCoords(blockId,coord)

        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

        ! Point to blocks center and face vars:
        call Grid_getBlkPtr(blockID,solnData,CENTER)
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)
        call Grid_getBlkPtr(blockID,facezData,FACEZ)

        ix1 = blkLimits(LOW,IAXIS)
        ix2 = blkLimits(HIGH,IAXIS) 

        jy1 = blkLimits(LOW,JAXIS)
        jy2 = blkLimits(HIGH,JAXIS)

        kz1 = blkLimits(LOW,KAXIS)
        kz2 = blkLimits(HIGH,KAXIS)


        facexData(VELC_FACE_VAR,ix1:ix2+1,jy1:jy2,kz1:kz2) = facexData(VELC_FACE_VAR,ix1:ix2+1,jy1:jy2,kz1:kz2) + &

                                                             (solnData(MDOT_VAR,ix1-1:ix2,jy1:jy2,kz1:kz2)+&
                                                              solnData(MDOT_VAR,ix1:ix2+1,jy1:jy2,kz1:kz2))/2.d0 * &

                                                             (solnData(NRMX_VAR,ix1-1:ix2,jy1:jy2,kz1:kz2)&
                                                              +solnData(NRMX_VAR,ix1:ix2+1,jy1:jy2,kz1:kz2))/2.d0 * &

                                                             ((solnData(ROLD_VAR,ix1-1:ix2,jy1:jy2,kz1:kz2)&
                                                              +solnData(ROLD_VAR,ix1:ix2+1,jy1:jy2,kz1:kz2))/2.d0 - &

                                                             (solnData(SMRH_VAR,ix1-1:ix2,jy1:jy2,kz1:kz2)&
                                                              +solnData(SMRH_VAR,ix1:ix2+1,jy1:jy2,kz1:kz2))/2.d0)


        faceyData(VELC_FACE_VAR,ix1:ix2,jy1:jy2+1,kz1:kz2) = faceyData(VELC_FACE_VAR,ix1:ix2,jy1:jy2+1,kz1:kz2) + &

                                                             (solnData(MDOT_VAR,ix1:ix2,jy1-1:jy2,kz1:kz2)+&
                                                              solnData(MDOT_VAR,ix1:ix2,jy1:jy2+1,kz1:kz2))/2.d0 * &

                                                             (solnData(NRMY_VAR,ix1:ix2,jy1-1:jy2,kz1:kz2)+&
                                                              solnData(NRMY_VAR,ix1:ix2,jy1:jy2+1,kz1:kz2))/2.d0 * &

                                                             ((solnData(ROLD_VAR,ix1:ix2,jy1-1:jy2,kz1:kz2)+&
                                                              solnData(ROLD_VAR,ix1:ix2,jy1:jy2+1,kz1:kz2))/2.d0 - &

                                                             (solnData(SMRH_VAR,ix1:ix2,jy1-1:jy2,kz1:kz2)+&
                                                              solnData(SMRH_VAR,ix1:ix2,jy1:jy2+1,kz1:kz2))/2.d0)

#if NDIM==3

        facezData(VELC_FACE_VAR,ix1:ix2,jy1:jy2,kz1:kz2+1) = facezData(VELC_FACE_VAR,ix1:ix2,jy1:jy2,kz1:kz2+1) + &
  
                                                             (solnData(MDOT_VAR,ix1:ix2,jy1:jy2,kz1-1:kz2)+&
                                                              solnData(MDOT_VAR,ix1:ix2,jy1:jy2,kz1:kz2+1))/2.d0 * &

                                                             (solnData(NRMZ_VAR,ix1:ix2,jy1:jy2,kz1-1:kz2)+&
                                                              solnData(NRMZ_VAR,ix1:ix2,jy1:jy2,kz1:kz2+1))/2.d0 * &

                                                             ((solnData(ROLD_VAR,ix1:ix2,jy1:jy2,kz1-1:kz2)+&
                                                              solnData(ROLD_VAR,ix1:ix2,jy1:jy2,kz1:kz2+1))/2.d0 - &

                                                             (solnData(SMRH_VAR,ix1:ix2,jy1:jy2,kz1-1:kz2)+&
                                                              solnData(SMRH_VAR,ix1:ix2,jy1:jy2,kz1:kz2+1))/2.d0)

#endif
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
  ins_predcorrflg = .true.
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)

  do lb = 1,blockCount

      blockID = blockList(lb)


        ! Get blocks dx, dy ,dz:
        call Grid_getDeltas(blockID,del)

        ! Get blocks coord and bsize
        ! Bounding box:
        call Grid_getBlkBoundBox(blockId,boundBox)
        bsize(1:NDIM) = boundBox(2,1:NDIM) - boundBox(1,1:NDIM)

        call Grid_getBlkCenterCoords(blockId,coord)

        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

        ! Point to blocks center and face vars:
        call Grid_getBlkPtr(blockID,solnData,CENTER)
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)
        call Grid_getBlkPtr(blockID,facezData,FACEZ)

        call ins_divergence_PC(facexData(VELC_FACE_VAR,:,:,:),&
                         faceyData(VELC_FACE_VAR,:,:,:),&
                         facezData(VELC_FACE_VAR,:,:,:),&
             blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
             blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
             blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
             blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS),&
             blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS),&
             blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS),&
                       del(DIR_X),del(DIR_Y),del(DIR_Z),&
                       solnData(DUST_VAR,:,:,:),&
                       solnData(DFUN_VAR,:,:,:),&
                       solnData(PFUN_VAR,:,:,:),&
                       solnData(NRMX_VAR,:,:,:),&
                       solnData(NRMY_VAR,:,:,:),&
                       solnData(NRMZ_VAR,:,:,:),&
                       solnData(SMRH_VAR,:,:,:),&
                       solnData(MDOT_VAR,:,:,:),&
                       mph_rho1,mph_rho2,&
                       facexData(RH1F_FACE_VAR,:,:,:),facexData(RH2F_FACE_VAR,:,:,:),&
                       faceyData(RH1F_FACE_VAR,:,:,:),faceyData(RH2F_FACE_VAR,:,:,:),&
                       facezData(RH1F_FACE_VAR,:,:,:),facezData(RH2F_FACE_VAR,:,:,:))


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

        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)


  end do

  poisfact = 1.0
  call Grid_solvePoisson (DELP_VAR, DUST_VAR, bc_types, bc_values, poisfact)

 do lb = 1,blockCount

      blockID = blockList(lb)


        ! Get blocks dx, dy ,dz:
        call Grid_getDeltas(blockID,del)

        ! Get blocks coord and bsize
        ! Bounding box:
        call Grid_getBlkBoundBox(blockId,boundBox)
        bsize(1:NDIM) = boundBox(2,1:NDIM) - boundBox(1,1:NDIM)

        call Grid_getBlkCenterCoords(blockId,coord)

        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

        ! Point to blocks center and face vars:
        call Grid_getBlkPtr(blockID,solnData,CENTER)
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)
        call Grid_getBlkPtr(blockID,facezData,FACEZ)

        sx = NGUARD+1
        sy = NGUARD*K2D+1
        sz = NGUARD*K3D+1
        ex = dataSize(DIR_X)-NGUARD
        ey = dataSize(DIR_Y)-NGUARD*K2D
        ez = dataSize(DIR_Z)-NGUARD*K3D

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

        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)


  end do

  gcMask = .FALSE.
  gcMask(NUNK_VARS+VELC_FACE_VAR) = .TRUE.                 ! u
  gcMask(NUNK_VARS+1*NFACE_VARS+VELC_FACE_VAR) = .TRUE.    ! v
#if NDIM == 3
  gcMask(NUNK_VARS+2*NFACE_VARS+VELC_FACE_VAR) = .TRUE.    ! w
#endif
  ins_predcorrflg = .false.
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)


end subroutine ins_MFcorrection
