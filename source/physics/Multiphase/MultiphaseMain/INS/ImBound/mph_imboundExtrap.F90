subroutine mph_imboundExtrap(blockCount, blockList,timeEndAdv,dt,dtOld,sweepOrder)

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

  real    :: hnorm2

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

  real :: dphidn

  real :: lambda_old(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD)
  real :: nx_mins, ny_mins, nz_mins, nx_plus, ny_plus, nz_plus
  real :: lambdax_mins, lambday_mins, lambdaz_mins, lambdax_plus, lambday_plus, lambdaz_plus
  real :: dt_ext
  integer :: extrap_index
  real :: nxconv, nyconv, nzconv

  do extrap_index=1,3
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

        lambda_old = solnData(LMDA_VAR,:,:,:)

        dt_ext = 0.5d0*del(IAXIS)

        k = 1
#if NDIM == 3
       do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
#endif
        do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
         do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
                         
            nxconv = -solnData(NMLX_VAR,i,j,k)
            nyconv = -solnData(NMLY_VAR,i,j,k)
            nzconv = -solnData(NMLZ_VAR,i,j,k)

            nx_plus = max(nxconv,0.)
            nx_mins = min(nxconv,0.)

            ny_plus = max(nyconv,0.)
            ny_mins = min(nyconv,0.)

            lambdax_plus = (lambda_old(i+1,j,k)-lambda_old(i,j,k))/del(IAXIS)
            lambday_plus = (lambda_old(i,j+1,k)-lambda_old(i,j,k))/del(JAXIS)

            lambdax_mins = (lambda_old(i,j,k)-lambda_old(i-1,j,k))/del(IAXIS)
            lambday_mins = (lambda_old(i,j,k)-lambda_old(i,j-1,k))/del(JAXIS)

#if NDIM == 3        
            nz_plus = max(nzconv,0.)
            nz_mins = min(nzconv,0.)

            lambdaz_plus = (lambda_old(i,j,k+1)-lambda_old(i,j,k))/del(IAXIS)
            lambdaz_mins = (lambda_old(i,j,k)-lambda_old(i,j,k-1))/del(JAXIS)
#else
            nz_plus = 0.0
            nz_mins = 0.0

            lambdaz_plus = 0.0
            lambdaz_mins = 0.0
#endif
           if(solnData(LMDA_VAR,i,j,k) .ge. 0.0) then

              solnData(LMDA_VAR,i,j,k) = lambda_old(i,j,k) - dt_ext*(nx_plus*lambdax_mins+nx_mins*lambdax_plus + &
                                                                     ny_plus*lambday_mins+ny_mins*lambday_plus + &
                                                                     nz_plus*lambdaz_mins+nz_mins*lambdaz_plus)

           end if
          
         end do
        end do
#if NDIM == 3
        end do
#endif

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

  end do
end subroutine mph_imboundExtrap
