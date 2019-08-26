subroutine mph_imbound(blockCount, blockList,timeEndAdv,dt,dtOld,sweepOrder)

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

  use Grid_interface, ONLY : Grid_getListOfBlocks, &
                             Grid_getDeltas,         &
                             Grid_getBlkBC,          &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkIndexLimits, &
                             Grid_fillGuardCells,    &
                             Grid_getBlkBoundBox,Grid_getBlkCenterCoords

  use IncompNS_data, ONLY : ins_alfa,ins_gravX,ins_gravY,ins_invRe,ins_gravZ

  use Multiphase_data, only: mph_rho1,mph_rho2,mph_sten,mph_crmx,mph_crmn, &
                             mph_vis1,mph_vis2,mph_lsit, mph_inls, mph_meshMe

  use Timers_interface, ONLY : Timers_start, Timers_stop

  use Driver_data, ONLY : dr_nstep, dr_simTime

  use ib_interface, ONLY : ib_stencils

  use ImBound_data, only : ib_stencil

  implicit none

#include "constants.h"
#include "IncompNS.h"
#include "ImBound.h"

  include "Flash_mpi.h"

  ! Arugments List
  integer, intent(in) :: sweepOrder
  integer, INTENT(INOUT) :: blockCount
  integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList
  real,    INTENT(IN) :: timeEndAdv,dt,dtOld

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  real, dimension(2,MDIM) :: boundBox


  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS), isAttached

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

  integer :: lb,blockID,ii,jj,kk,ierr,i,j,k,dir

  real bsize(MDIM),coord(MDIM), vel_probe(MDIM)
  
  real del(MDIM),xcell,ycell,zcell,rc

  real    :: r_avg
  integer :: n_avg
  !kpd
  real :: lsDT,lsT,minCellDiag
  real :: volSum,volSumAll

  real :: vol, cx, cy, vx, vy
  real :: xh, yh, xl, yl

  !- kpd - For Overall Solver Timer... 
  real :: t_startMP1,t_stopMP1,t_startMP2,t_startMP2a,t_stopMP2

  !- kpd - For Poisson Timer...
  integer :: t

  integer :: listofBlocks(MAXBLOCKS)
  integer :: count
  integer :: intval

  real    :: hnorm, xprobe(3), yprobe(3), zprobe(3), phiprobe

  real,parameter  :: htol = 0.0001

  integer :: gridfl(MDIM)
  real    :: externalPt(MDIM), part_Nml(MDIM), dfe
  integer, dimension(ib_stencil,MDIM) :: ib_external
  real, dimension(ib_stencil,NDIM+1) :: ib_external_phile
  integer,parameter,dimension(MDIM):: FACE_IND =(/FACEX,FACEY,FACEZ/)

  integer :: idim
  real    :: xyz_stencil(ib_stencil,MDIM)
  real :: delaux(MDIM)
  real :: zp(3)

  integer, parameter :: derivflag = 0
  integer :: ib_ind

  real    :: hratio, veli, this_psi

  real    :: nrmx, nrmy, nmlx, nmly, ib_theta

  integer :: probe_index

  do lb = 1,blockCount

        blockID = blockList(lb)

        call Grid_getBlkBoundBox(blockId,boundBox)

        bsize(:) = boundBox(2,:) - boundBox(1,:)

        call Grid_getBlkCenterCoords(blockId,coord)
 
        call Grid_getDeltas(blockID,del)

        ! Get Blocks internal limits indexes:
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

        ! Point to blocks center and face vars:
        call Grid_getBlkPtr(blockID,solnData,CENTER)
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)
        call Grid_getBlkPtr(blockID,facezData,FACEZ)

        k = 1
#if NDIM == 3
       do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
#endif
        do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
         do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
                         
           xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
                   real(i - NGUARD - 1)*del(IAXIS) +   &
                   0.5*del(IAXIS)

           ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                   real(j - NGUARD - 1)*del(JAXIS)  +  &
                   0.5*del(JAXIS)

          zcell = 0.0

#if NDIM == 3
           zcell  = coord(KAXIS) - bsize(KAXIS)/2.0 +  &
                   real(k - NGUARD - 1)*del(KAXIS)  +  &
                   0.5*del(KAXIS)
#endif
          
           if(solnData(LMDA_VAR,i,j,k) .ge. 0.0 .and. solnData(LMDA_VAR,i,j,k) .le. 1.5*del(IAXIS)) then

           ! Get probe in fluid
           hnorm = 1.0*del(JAXIS)

           xprobe(1) = xcell + solnData(NMLX_VAR,i,j,k)*(solnData(LMDA_VAR,i,j,k)+hnorm)
           yprobe(1) = ycell + solnData(NMLY_VAR,i,j,k)*(solnData(LMDA_VAR,i,j,k)+hnorm)
           zprobe(1) = 0.0

#if NDIM == 3
           zprobe(1) = zcell + solnData(NMLZ_VAR,i,j,k)*(solnData(LMDA_VAR,i,j,k)+hnorm)
#endif
           !xprobe(2) = xprobe(1) + solnData(TNGX_VAR,i,j,k)*del(IAXIS)
           !yprobe(2) = yprobe(1) + solnData(TNGY_VAR,i,j,k)*del(JAXIS)
           !zprobe(2) = 0.0

           !xprobe(3) = xprobe(1) - solnData(TNGX_VAR,i,j,k)*del(IAXIS)
           !yprobe(3) = yprobe(1) - solnData(TNGY_VAR,i,j,k)*del(JAXIS)
           !zprobe(3) = 0.0

           ! Interpolate function at probe 
           do probe_index = 1,1
           externalPt(IAXIS) = xprobe(probe_index)
           externalPt(JAXIS) = yprobe(probe_index)
           externalPt(KAXIS) = zprobe(probe_index)

           part_Nml(IAXIS) = solnData(NMLX_VAR,i,j,k)
           part_Nml(JAXIS) = solnData(NMLY_VAR,i,j,k)
           part_Nml(KAXIS) = 0.0

#if NDIM == 3
           part_Nml(KAXIS) = solnData(NMLZ_VAR,i,j,k)
#endif

           ! Cell centered stencil for DFUN interpolation at probe
           gridfl(:) = CENTER

           call ib_stencils(externalPt,part_Nml,gridfl,del,coord,bsize, &
                            ib_external(:,:),dfe,FORCE_FLOW)

           delaux = 0.5*del

           xyz_stencil(:,:) = 0.

           do idim = 1,NDIM
              xyz_stencil(:,idim) = coord(idim) - 0.5*bsize(idim) + &
                                    real(ib_external(1:ib_stencil,idim) - NGUARD - 1)*del(idim) + delaux(idim)
           enddo

           ! Get shape function ib_external_phile  for points on stencil
           call ib_getInterpFunc(externalPt,xyz_stencil,del,derivflag,ib_external_phile)

           zp(probe_index) = 0.      ! zp = DFUN at probe point

#if NDIM == 2
           do ib_ind = 1 , ib_stencil
                zp(probe_index) = zp(probe_index) + ib_external_phile(ib_ind,CONSTANT_ONE) * &
                solnData(DFUN_VAR,ib_external(ib_ind,IAXIS),ib_external(ib_ind,JAXIS),1);
           enddo
#endif

#if NDIM == 3
           do ib_ind = 1 , ib_stencil
                zp(probe_index) = zp(probe_index) + ib_external_phile(ib_ind,CONSTANT_ONE) * &
                solnData(DFUN_VAR,ib_external(ib_ind,IAXIS),ib_external(ib_ind,JAXIS),ib_external(ib_ind,KAXIS));
           enddo
#endif

           enddo

           hratio = (solnData(LMDA_VAR,i,j,k) + hnorm)

           !if(zp(1)*zp(2) .le. 0.0 .or. zp(1)*zp(3) .le. 0.0 .or. zp(1) .ge. 0.0) then
           solnData(DFUN_VAR,i,j,k) = zp(1)-hratio*cos(60.0*acos(-1.0)/180)
           !end if

           end if
          
         end do
        end do
#if NDIM == 3
        end do
#endif

!        k = 1
!#if NDIM == 3
!       do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
!#endif
!        do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
!         do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
! 
!          if(solnData(LMDA_VAR,i,j,k) .gt. 1.5*del(IAXIS)) &
!             solnData(DFUN_VAR,i,j,k) = &
!             min(-solnData(LMDA_VAR,i,j,k) + 1.5*del(IAXIS),solnData(DFUN_VAR,i,j,k))
!
!         end do
!        end do
!#if NDIM == 3
!        end do
!#endif
        ! Release pointers:
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  end do
 
  gcMask = .FALSE.

  gcMask(DFUN_VAR) = .TRUE.

  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)

end subroutine mph_imbound
