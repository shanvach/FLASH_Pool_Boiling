subroutine Heat_getIBFlux(blockCount, blockList,timeEndAdv,dt,dtOld,sweepOrder)

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
                             mph_vis1,mph_vis2,mph_lsit, mph_inls, mph_meshMe, &
                             mph_thco1, mph_thco2

  use Timers_interface, ONLY : Timers_start, Timers_stop

  use Driver_data, ONLY : dr_nstep, dr_simTime

  use ib_interface, ONLY : ib_stencils

  use ImBound_data, only : ib_stencil

  use Heat_AD_interface, only: Heat_GFMstencil_o1

  use Heat_AD_data, only: ht_hflux_flag, ht_hflux_counter, &
                          ht_ibx, ht_iby, ht_ibz, ht_ibT, ht_ibNu, &
                          ht_ibhflux_counter, ht_ibNu_t, ht_ibNu_l
 
  use Simulation_data, only: sim_xMin, sim_xMax, sim_yMin, sim_yMax
 
  implicit none

#include "constants.h"
#include "Heat_AD.h"
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

  real    :: hnorm(3), xprobe(3), yprobe(3), zprobe(3), phiprobe(3)

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

  real :: thxp, thxm, thyp, thym, thzp, thzm

  real :: Tboundary = 1.0
 
  real :: tol = 0.1
 
  real :: thco3, thcoE, thcoP, thcoC

  integer :: hflux_counter

  real :: h2,h3, hflux_buffer

hflux_buffer = 1.0

#define AROUND_INTERFACE 1

  thco3 = 10*mph_thco2 

   if(ht_hflux_flag) then

    ht_hflux_counter = 0

    do lb = 1,blockCount

     blockID = blockList(lb)
     call Grid_getDeltas(blockID,del)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkBoundBox(blockId,boundBox)
     bsize(:) = boundBox(2,:) - boundBox(1,:)
     call Grid_getBlkCenterCoords(blockId,coord)

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

#if NDIM == 2
#if AROUND_INTERFACE == 1
           if((solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i+1,j,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i-1,j,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i,j+1,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i,j-1,k) .le. 0.0 ) .and. &
               solnData(LMDA_VAR,i,j,k) .ge. 0.0 .and. abs(solnData(DFUN_VAR,i,j,k)) .le. hflux_buffer) &
               ht_hflux_counter = ht_hflux_counter + 1
#else
           if((solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i+1,j,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i-1,j,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i,j+1,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i,j-1,k) .le. 0.0 ) .and. &
               solnData(LMDA_VAR,i,j,k) .ge. 0.0) ht_hflux_counter = ht_hflux_counter + 1
#endif

#else
#if AROUND_INTERFACE == 1
           if((solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i+1,j,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i-1,j,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i,j+1,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i,j-1,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i,j,k+1) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i,j,k-1) .le. 0.0 ) .and. & 
               solnData(LMDA_VAR,i,j,k) .ge. 0.0 .and. abs(solnData(DFUN_VAR,i,j,k)) .le. hflux_buffer) &
               ht_hflux_counter = ht_hflux_counter + 1
#else
           if((solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i+1,j,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i-1,j,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i,j+1,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i,j-1,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i,j,k+1) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i,j,k-1) .le. 0.0 ) .and. & 
               solnData(LMDA_VAR,i,j,k) .ge. 0.0) ht_hflux_counter = ht_hflux_counter + 1
#endif
#endif
           end do
        end do
#if NDIM == 3
        end do
#endif

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)

    end do

    allocate(ht_ibx(ht_hflux_counter), ht_iby(ht_hflux_counter), ht_ibz(ht_hflux_counter))
    allocate(ht_ibT(ht_hflux_counter), ht_ibNu(ht_hflux_counter))
    
    ht_hflux_flag = .false.

  end if

  hflux_counter = 1

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
         
#if NDIM == 2
#if AROUND_INTERFACE == 1
           if((solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i+1,j,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i-1,j,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i,j+1,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i,j-1,k) .le. 0.0 ) .and. &
               solnData(LMDA_VAR,i,j,k) .ge. 0.0 .and. abs(solnData(DFUN_VAR,i,j,k)) .le. hflux_buffer) then
#else
           if((solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i+1,j,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i-1,j,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i,j+1,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i,j-1,k) .le. 0.0 ) .and. &
               solnData(LMDA_VAR,i,j,k) .ge. 0.0) then
#endif

#else
#if AROUND_INTERFACE == 1
           if((solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i+1,j,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i-1,j,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i,j+1,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i,j-1,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i,j,k+1) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i,j,k-1) .le. 0.0 ) .and. & 
               solnData(LMDA_VAR,i,j,k) .ge. 0.0 .and. abs(solnData(DFUN_VAR,i,j,k)) .le. hflux_buffer) then
#else
           if((solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i+1,j,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i-1,j,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i,j+1,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i,j-1,k) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i,j,k+1) .le. 0.0 .or.  &
               solnData(LMDA_VAR,i,j,k)*solnData(LMDA_VAR,i,j,k-1) .le. 0.0 ) .and. & 
               solnData(LMDA_VAR,i,j,k) .ge. 0.0) then
#endif
#endif

           ! Get probe in fluid
           hnorm(1) = 0.0
           hnorm(2) = sign(1.5*del(IAXIS),solnData(LMDA_VAR,i,j,k))
           hnorm(3) = sign(2.0*del(IAXIS),solnData(LMDA_VAR,i,j,k))
               
           do probe_index = 1,3

           xprobe(probe_index) = xcell + solnData(NMLX_VAR,i,j,k)*(solnData(LMDA_VAR,i,j,k)+hnorm(probe_index))
           yprobe(probe_index) = ycell + solnData(NMLY_VAR,i,j,k)*(solnData(LMDA_VAR,i,j,k)+hnorm(probe_index))
           zprobe(probe_index) = 0.0

#if NDIM == 3
           zprobe(probe_index) = zcell + solnData(NMLZ_VAR,i,j,k)*(solnData(LMDA_VAR,i,j,k)+hnorm(probe_index))
#endif

           end do

           ! Interpolate function at probe 
           do probe_index = 1,3
           externalPt(IAXIS) = xprobe(probe_index)
           externalPt(JAXIS) = yprobe(probe_index)
           externalPt(KAXIS) = zprobe(probe_index)

           part_Nml(IAXIS) = solnData(NMLX_VAR,i,j,k)
           part_Nml(JAXIS) = solnData(NMLY_VAR,i,j,k)
           part_Nml(KAXIS) = 0.0

#if NDIM == 3
           part_Nml(KAXIS) = solnData(NMLZ_VAR,i,j,k)
#endif

           ! Cell centered stencil for temperature interpolation at probe
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

           zp(probe_index)  = 0.      ! zp = temperature at probe point
           phiprobe(probe_index) = 0.
#if NDIM == 2
           do ib_ind = 1 , ib_stencil
                zp(probe_index) = zp(probe_index) + ib_external_phile(ib_ind,CONSTANT_ONE) * &
                solnData(TEMP_VAR,ib_external(ib_ind,IAXIS),ib_external(ib_ind,JAXIS),1);

                phiprobe(probe_index) = phiprobe(probe_index) + ib_external_phile(ib_ind,CONSTANT_ONE) * &
                solnData(DFUN_VAR,ib_external(ib_ind,IAXIS),ib_external(ib_ind,JAXIS),1);
           enddo
#endif

#if NDIM == 3
           do ib_ind = 1 , ib_stencil
                zp(probe_index) = zp(probe_index) + ib_external_phile(ib_ind,CONSTANT_ONE) * &
                solnData(TEMP_VAR,ib_external(ib_ind,IAXIS),ib_external(ib_ind,JAXIS),ib_external(ib_ind,KAXIS));
           enddo

           do ib_ind = 1 , ib_stencil
                phiprobe(probe_index) = phiprobe(probe_index) + ib_external_phile(ib_ind,CONSTANT_ONE) * &
                solnData(DFUN_VAR,ib_external(ib_ind,IAXIS),ib_external(ib_ind,JAXIS),ib_external(ib_ind,KAXIS));
           enddo
#endif

           enddo
        
           if(solnData(LMDA_VAR,i,j,k) .lt. 0.0) then

                thcoC = mph_thco1/mph_thco2
                if(solnData(DFUN_VAR,i,j,k) .lt. 0.0) thcoC = mph_thco2/mph_thco2

                thcoP = thco3/mph_thco2
        
           else

                thcoC = thco3/mph_thco2

                thcoP = mph_thco1/mph_thco2
                if(phiprobe(2) .lt. 0.0) thcoP = mph_thco2/mph_thco2


           end if
        
           thcoE = ((solnData(LMDA_VAR,i,j,k)+hnorm(2))*thcoC*thcoP)/(thcoP*solnData(LMDA_VAR,i,j,k) + thcoC*hnorm(2))

           hratio = thcoE/thcoP

           ht_ibx(hflux_counter)  = xprobe(1)
           ht_iby(hflux_counter)  = yprobe(1)
           ht_ibz(hflux_counter)  = zprobe(1)

           ht_ibT(hflux_counter)  = zp(1)

           ! First Order Approximation
           ht_ibNu(hflux_counter) = hratio*(zp(2)-solnData(TEMP_VAR,i,j,k))/(solnData(LMDA_VAR,i,j,k)+hnorm(2))


           ! Second Order Approximation
           !h2 = solnData(LMDA_VAR,i,j,k) + hnorm(2)
           !h3 = solnData(LMDA_VAR,i,j,k) + hnorm(3)

           !ht_ibNu(hflux_counter) = hratio*(zp(3)*h2*h2 - zp(2)*h3*h3 - &
           !                                 solnData(TEMP_VAR,i,j,k)*(h2*h2 - h3*h3))/(h2*h3*(h2-h3))          

           hflux_counter = hflux_counter+1

           end if
           !---------------------------------------------------------------------!

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

  ht_ibNu_t = sum(ht_ibNu)

end subroutine Heat_getIBFlux
