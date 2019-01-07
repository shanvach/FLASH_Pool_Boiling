subroutine mph_imboundJumps(blockCount, blockList,timeEndAdv,dt,dtOld,sweepOrder)

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
                             mph_vis1,mph_vis2,mph_lsit, mph_inls, mph_meshMe,&
                             mph_radius, mph_isAttached, mph_timeStamp, mph_vlim, mph_psi_adv

  use Timers_interface, ONLY : Timers_start, Timers_stop

  use Driver_data, ONLY : dr_nstep, dr_simTime

  use ib_interface, ONLY : ib_stencils

  use ImBound_data, only : ib_stencil

  use Heat_AD_data, only : ht_psi

  ! Following routine is written by Akash
  ! Actual calls written by Shizao and Keegan
  ! This subroutine decouples Multiphase calls from ins_ab2rk3_VD 

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

  real    :: hnorm, xprobe, yprobe, zprobe, phiprobe

  real,parameter  :: htol = 0.0001

  integer :: gridfl(MDIM)
  real    :: externalPt(MDIM), part_Nml(MDIM), dfe
  integer, dimension(ib_stencil,MDIM) :: ib_external
  real, dimension(ib_stencil,NDIM+1) :: ib_external_phile
  integer,parameter,dimension(MDIM):: FACE_IND =(/FACEX,FACEY,FACEZ/)

  integer :: idim
  real    :: xyz_stencil(ib_stencil,MDIM)
  real :: delaux(MDIM)
  real :: zp, lambdaX, lambdaY

  integer, parameter :: derivflag = 0
  integer :: ib_ind

  real    :: hratio, veli, this_psi

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
        do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
         do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

           if(solnData(LMDA_VAR,i,j,k) .ge. 0.0 .and. solnData(LMDA_VAR,i,j,k) .lt. 2.0*del(IAXIS)) then
                          
           xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
                   real(i - NGUARD - 1)*del(IAXIS) +   &
                   0.5*del(IAXIS)

           ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                   real(j - NGUARD - 1)*del(JAXIS)  +  &
                   0.5*del(JAXIS)

           !zcell  = coord(KAXIS) - bsize(KAXIS)/2.0 +  &
           !        real(k - NGUARD - 1)*del(KAXIS)  +  &
           !        0.5*del(KAXIS)
          
           ! Get probe in fluid
           hnorm   = 1.*sqrt((1.2*del(IAXIS)*solnData(NMLX_VAR,i,j,k))**2. + (1.2*del(JAXIS)*solnData(NMLY_VAR,i,j,k))**2.)

           xprobe = xcell + solnData(NMLX_VAR,i,j,k)*(solnData(LMDA_VAR,i,j,k)+hnorm)
           yprobe = ycell + solnData(NMLY_VAR,i,j,k)*(solnData(LMDA_VAR,i,j,k)+hnorm)

           ! Interpolate function at probe 
           externalPt(IAXIS) = xprobe
           externalPt(JAXIS) = yprobe
           externalPt(KAXIS) = 0.0

           part_Nml(IAXIS) = solnData(NMLX_VAR,i,j,k)
           part_Nml(JAXIS) = solnData(NMLY_VAR,i,j,k)
           part_Nml(KAXIS) = 0.0

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

           zp = 0.      ! zp = DFUN at probe point

           do ib_ind = 1 , ib_stencil
                zp = zp + ib_external_phile(ib_ind,CONSTANT_ONE) * &
                solnData(SIGP_VAR,ib_external(ib_ind,IAXIS),ib_external(ib_ind,JAXIS),1);
           enddo

           do dir=1,NDIM

                gridfl(:) = CENTER
                gridfl(dir) = FACES
                delaux(1:NDIM)   = 0.5*del(1:NDIM)
                delaux(dir) = 0.

                ! Define Interpolation Stencil For Particle:
                call ib_stencils(externalPt,part_Nml,gridfl,del,coord,bsize,   &
                                 ib_external(:,:),dfe,FORCE_FLOW)

                ! Interpolation of the values of velocity to Lagrangian points:
                xyz_stencil(:,:) = 0. 
                do idim = 1,NDIM
                    xyz_stencil(:,idim) = coord(idim) - 0.5*bsize(idim) + &
                        real(ib_external(1:ib_stencil,idim) - NGUARD - 1)*del(idim) + delaux(idim) 
                enddo
                call ib_getInterpFunc(externalPt,xyz_stencil,del,0,ib_external_phile)
                vel_probe(dir) = 0.
                do ii = 1 , ib_stencil      
                    select case(dir)
                    case(FACEX) 
                    vel_probe(dir) = vel_probe(dir) + ib_external_phile(ii,CONSTANT_ONE) * &
                            facexData(SIGM_FACE_VAR,ib_external(ii,IAXIS),ib_external(ii,JAXIS),ib_external(ii,KAXIS));   
                    case(FACEY) 
                    vel_probe(dir) = vel_probe(dir) + ib_external_phile(ii,CONSTANT_ONE) * &
                            faceyData(SIGM_FACE_VAR,ib_external(ii,IAXIS),ib_external(ii,JAXIS),ib_external(ii,KAXIS));   
#if NDIM == MDIM
                    case(FACEZ)
                    vel_probe(dir) = vel_probe(dir) + ib_external_phile(ii,CONSTANT_ONE) * &
                            facezData(SIGM_FACE_VAR,ib_external(ii,IAXIS),ib_external(ii,JAXIS),ib_external(ii,KAXIS));   
#endif
                    end select
                enddo
            end do

             solnData(SIGP_VAR,i,j,k) = zp

             lambdaX = 0.5*(solnData(LMDA_VAR,i,j,k)+solnData(LMDA_VAR,i-1,j,k))
             lambdaY = 0.5*(solnData(LMDA_VAR,i,j,k)+solnData(LMDA_VAR,i,j-1,k))

             facexData(SIGM_FACE_VAR,i,j,k) = vel_probe(1)
             faceyData(SIGM_FACE_VAR,i,j,k) = vel_probe(2)

           end if 

         end do
        end do

         ! Release pointers:
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  end do
 
  !gcMask = .TRUE.

  !call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
  !     maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)

end subroutine mph_imboundJumps
