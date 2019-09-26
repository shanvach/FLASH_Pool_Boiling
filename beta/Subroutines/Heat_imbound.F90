subroutine Heat_imbound(blockCount,blockList,timeEndAdv,dt,ivar)

#include "Flash.h"

!#define GCELL_FORCING
#define LINE_FORCING

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

  use IncompNS_data, ONLY : ins_alfa,ins_gravX,ins_gravY,ins_invRe,ins_gravZ,&
                            ins_predcorrflg

  use Multiphase_data, only: mph_rho1,mph_rho2,mph_sten,mph_crmx,mph_crmn, &
                             mph_vis1,mph_vis2,mph_lsit, mph_inls, mph_meshMe,&
                             mph_radius, mph_isAttached, mph_timeStamp, mph_vlim, mph_psi_adv

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
  integer, INTENT(INOUT) :: blockCount
  integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList
  real,    INTENT(IN) :: timeEndAdv,dt
  integer, intent(in) :: ivar

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

  real :: vol, cx, cy
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
  real :: zp

  integer, parameter :: derivflag = 0
  integer :: ib_ind

  real    :: hratio, temp, this_psi

  real :: lambda
  real :: nx, ny, m, n

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

           lambda = solnData(LMDA_VAR,i,j,k)

           nx = solnData(NMLX_VAR,i,j,k)
           ny = solnData(NMLY_VAR,i,j,k)

#ifdef LINE_FORCING
           if(abs(lambda) .le. 1.0*del(IAXIS)) then
#endif

#ifdef GCELL_FORCING
           if(lambda .ge. 0.0 .and. lambda .le. 1.0*del(IAXIS)) then
#endif
               
           xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
                   real(i - NGUARD - 1)*del(IAXIS) +   &
                   0.5*del(IAXIS)

           ycell = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                   real(j - NGUARD - 1)*del(JAXIS) +  &
                   0.5*del(JAXIS)
         
           ! Get probe in fluid

#ifdef LINE_FORCING
           if(lambda .lt. 0.0) then
             hnorm =  2.0*del(IAXIS)
           else
             hnorm = -2.0*del(IAXIS)
           end if
#endif

#ifdef GCELL_FORCING
           hnorm = 1.0*del(IAXIS)
#endif

           xprobe = xcell + nx*(lambda+hnorm)
           yprobe = ycell + ny*(lambda+hnorm)

           ! Interpolate function at probe 
           externalPt(IAXIS) = xprobe
           externalPt(JAXIS) = yprobe
           externalPt(KAXIS) = 0.0

           part_Nml(IAXIS) = nx
           part_Nml(JAXIS) = ny
           part_Nml(KAXIS) = 0.0

           gridfl(:)   = CENTER
           delaux(1:NDIM) = 0.5*del(1:NDIM)

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
                
           temp = 0

           do ii = 1 , ib_stencil      
              temp = temp + ib_external_phile(ii,CONSTANT_ONE) * &
              solnData(ivar,ib_external(ii,IAXIS),ib_external(ii,JAXIS),ib_external(ii,KAXIS));   
           end do

           m = abs(lambda)
           n = abs(lambda+hnorm)

#ifdef LINE_FORCING
           solnData(ivar,i,j,k) = (m*temp+n)/(m+n)    
#endif

#ifdef GCELL_FORCING
           solnData(ivar,i,j,k) = (2.0*(lambda+hnorm) - temp*lambda)/hnorm
#endif

           end if

         end do
        end do

         ! Release pointers:
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  end do
 
  gcMask = .FALSE.
  gcMask(TEMP_VAR) = .TRUE.

  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)

end subroutine Heat_imbound
