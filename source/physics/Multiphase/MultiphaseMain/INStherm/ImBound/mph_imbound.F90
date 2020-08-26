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
                             mph_vis1,mph_vis2,mph_lsit, mph_inls, mph_meshMe,&
                             mph_vlim, mph_psi_adv

  use Timers_interface, ONLY : Timers_start, Timers_stop

  use Driver_data, ONLY : dr_nstep, dr_simTime

  use ib_interface, ONLY : ib_stencils

  use ImBound_data, only : ib_stencil

  use Heat_AD_data, only : ht_psi

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
  real    :: part_Tng(MDIM)
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
          
           if(solnData(LMDA_VAR,i,j,k) .gt. 0.0 .and. solnData(LMDA_VAR,i,j,k) .le. 1.5*del(IAXIS)) then

           ! Get probe in fluid
           hnorm  = 1.0*del(JAXIS)

           xprobe(1) = xcell + solnData(NMLX_VAR,i,j,k)*(solnData(LMDA_VAR,i,j,k)+hnorm)
           yprobe(1) = ycell + solnData(NMLY_VAR,i,j,k)*(solnData(LMDA_VAR,i,j,k)+hnorm)
           zprobe(1) = 0.0

           xprobe(2) = xprobe(1) + solnData(TNGX_VAR,i,j,k)*hnorm
           yprobe(2) = yprobe(1) + solnData(TNGY_VAR,i,j,k)*hnorm
           zprobe(2) = 0.0

           xprobe(3) = xprobe(1) - solnData(TNGX_VAR,i,j,k)*hnorm
           yprobe(3) = yprobe(1) - solnData(TNGY_VAR,i,j,k)*hnorm
           zprobe(3) = 0.0

#if NDIM == 3
           zprobe(1) = zcell + solnData(NMLZ_VAR,i,j,k)*(solnData(LMDA_VAR,i,j,k)+hnorm)
           zprobe(2) = zprobe(1) + solnData(TNGZ_VAR,i,j,k)*hnorm
           zprobe(3) = zprobe(1) - solnData(TNGZ_VAR,i,j,k)*hnorm
#endif
           ! Interpolate function at probe 
           do probe_index = 1,3
           externalPt(IAXIS) = xprobe(probe_index)
           externalPt(JAXIS) = yprobe(probe_index)
           externalPt(KAXIS) = zprobe(probe_index)

           part_Nml(IAXIS) = solnData(NMLX_VAR,i,j,k)
           part_Nml(JAXIS) = solnData(NMLY_VAR,i,j,k)
           part_Nml(KAXIS) = 0.0

           part_Tng(IAXIS) = solnData(NRMX_VAR,i,j,k)
           part_Tng(JAXIS) = solnData(NRMY_VAR,i,j,k)
           part_Tng(KAXIS) = 0.0

#if NDIM == 3
           part_Nml(KAXIS) = solnData(NMLZ_VAR,i,j,k)
           part_Tng(KAXIS) = solnData(NRMZ_VAR,i,j,k)
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


           if(probe_index == 1) then

           veli=0

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

                call ib_getInterpFunc(externalPt,xyz_stencil,del,derivflag,ib_external_phile)

                vel_probe(dir) = 0.

                do ib_ind = 1 , ib_stencil      
                    select case(dir)
                    case(FACEX) 
                    vel_probe(dir) = vel_probe(dir) + ib_external_phile(ib_ind,CONSTANT_ONE) * &
                            facexData(VELI_FACE_VAR,ib_external(ib_ind,IAXIS),ib_external(ib_ind,JAXIS),ib_external(ib_ind,KAXIS));   
                    case(FACEY) 
                    vel_probe(dir) = vel_probe(dir) + ib_external_phile(ib_ind,CONSTANT_ONE) * &
                            faceyData(VELI_FACE_VAR,ib_external(ib_ind,IAXIS),ib_external(ib_ind,JAXIS),ib_external(ib_ind,KAXIS));   
#if NDIM == MDIM
                     case(FACEZ)
                    vel_probe(dir) = vel_probe(dir) + ib_external_phile(ib_ind,CONSTANT_ONE) * &
                            facezData(VELI_FACE_VAR,ib_external(ib_ind,IAXIS),ib_external(ib_ind,JAXIS),ib_external(ib_ind,KAXIS));   
#endif
                    end select
                enddo 

               veli = veli + vel_probe(dir) * part_Tng(dir)

           enddo
           end if

           enddo

            this_psi = ht_psi

            !if(zp(1)*zp(2) .le. 0.0 .or. zp(1)*zp(3) .le. 0.0) then
            !if(veli .ge. 0.0) then
            !     if(abs(veli) .le. 0.2) then
            !          this_psi = ((mph_psi_adv - ht_psi)/(2*mph_vlim))*abs(veli)+ &
            !                                  (mph_psi_adv + ht_psi)/2.0d0
            !
            !     else
            !     this_psi = mph_psi_adv
            !    
            !     end if
            !end if
            !end if

           dphidn = cos(this_psi)
           !dphidn = cos(90*acos(-1.0)/180)
           !dphidn = (zp(2)-zp(1))/(hnorm2-hnorm)

           hratio = (solnData(LMDA_VAR,i,j,k) + hnorm)
           solnData(DFUN_VAR,i,j,k) = zp(1) - hratio*dphidn

           end if
          
         end do
        end do
#if NDIM == 3
        end do
#endif

        k = 1
#if NDIM == 3
       do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
#endif
        do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
         do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
 
          if(solnData(LMDA_VAR,i,j,k) .gt. 1.5*del(IAXIS)) &
             solnData(DFUN_VAR,i,j,k) = &
             min(-solnData(LMDA_VAR,i,j,k) + 1.5*del(IAXIS),solnData(DFUN_VAR,i,j,k))

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

end subroutine mph_imbound
