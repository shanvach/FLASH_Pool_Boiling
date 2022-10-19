subroutine mph_advect(blockCount, blockList, timeEndAdv, dt,dtOld,sweepOrder)

#define NUCLEATE_BOILING
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

  use IncompNS_data, ONLY : ins_meshMe,ins_alfa,ins_gravX,ins_gravY,ins_invRe,ins_gravZ

  use Multiphase_data, only: mph_rho1,mph_rho2,mph_sten,mph_crmx,mph_crmn, &
                             mph_vis1,mph_vis2,mph_lsit, mph_inls, mph_meshMe,&
                             mph_radius, mph_isAttached, mph_timeStamp, &
                             mph_isAttachedAll, mph_timeStampAll,&
                             mph_isAttachedOld, mph_nucSiteTemp, mph_offset, mph_psi, mph_psi_adv, mph_vlim

  use mph_interface, only : mph_KPDcurvature2DAB, mph_KPDcurvature2DC, &
                            mph_KPDadvectWENO3, mph_KPDlsRedistance,  &
                            mph_KPDcurvature3DAB, mph_KPDcurvature3DC,&
                            mph_KPDadvectWENO3_3D, mph_KPDlsRedistance_3D,&
                            mph_getInterfaceVelocity,mph_getInterfaceVelocity_3D 


  use Timers_interface, ONLY : Timers_start, Timers_stop

  use Driver_data, ONLY : dr_nstep, dr_simTime

!  use Heat_AD_data, ONLY: ht_psi, ht_tWait, ht_Tnuc

  use Simulation_data, ONLY: sim_nuc_site_x, sim_nuc_site_y, sim_nuc_radii, sim_nuc_site_z, sim_nucSiteDens
!------------------------------------------------------------------------------------------------------------------Shantanu
  use Heat_AD_data, ONLY: ht_psi, ht_tWait_hydrophilic, ht_tWait_hydrophobic, ht_Tnuc 

 use Simulation_data, ONLY : sim_theta_hydrophobic, sim_theta_hydrophilic,hydrophobicity_ratio,sim_x_start_point,sim_x_end_point,sim_z_start_point,sim_z_end_point

!-------------------------------------------------------------------------------------------------------------------
  ! Following routine is written by Akash
  ! Actual calls written by Shizao and Keegan
  ! This subroutine decouples Multiphase calls from ins_ab2rk3_VD 

  implicit none

#include "constants.h"
#include "IncompNS.h"

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

  integer :: lb,blockID,ii,jj,kk,ierr,i,j,k

  real bsize(MDIM),coord(MDIM)
  
  real del(MDIM),xcell,ycell,zcell,rc,xcellp,zcellp

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

  integer :: nuc_index, tSI
  real    :: nuc_dfun, nucSiteTemp

  real    :: tol=1E-13

  real    :: veli
!--------------------------------------------------------------------------------------------------Shantanu

 ! real    :: theta_hydrophobic = (165.0/180.0)*acos(-1.0)
 ! real    :: theta_hydrophilic = (20.0/180.0)*acos(-1.0)
 ! real    :: pitch  = 0.25
 ! real    :: dia    = 0.75
 ! real    :: start_point      = -8
  real    :: x_position,z_position
  integer :: counter,nspots,count_x,count_z
  real,dimension(:),allocatable :: x_left,x_right,z_bottom,z_top
 real, dimension(sim_nucSiteDens) :: waittime

  waittime(:) =  ht_tWait_hydrophilic

!  com_dif  = sim_pitch + sim_dia
  nspots    = int((sim_x_end_point - sim_x_start_point))
  allocate(x_left(nspots**2))
  allocate(x_right(nspots**2))
  allocate(z_top(nspots**2))
  allocate(z_bottom(nspots**2))

  x_left(1)  = sim_x_start_point + 0.5*(1-sqrt(hydrophobicity_ratio))
 ! x_right(1) = x_left(1) + sqrt(hydrophobicity_ratio)
  z_bottom(1)  = sim_z_start_point +  0.5*(1-sqrt(hydrophobicity_ratio))
 ! z_top(1) =  z_bottom(1) +  sqrt(hydrophobicity_ratio)
 do count_z = 1,nspots
        x_left((count_z -1)*nspots+1) = x_left(1)
        z_bottom((count_z -1)*nspots+1) = z_bottom(1) + (count_z-1)
        do count_x = 1,nspots
                x_left((count_z -1)*nspots + count_x)  =  x_left((count_z -1)*nspots+1) + (count_x-1)
                z_bottom((count_z -1)*nspots + count_x)  = z_bottom((count_z -1)*nspots+1)
        end do
  end do
  x_right(:)  = x_left(:) + sqrt(hydrophobicity_ratio)
  z_top(:)  = z_bottom(:) + sqrt(hydrophobicity_ratio)


            do nuc_index = 1,sim_nucSiteDens
                do counter = 1,nspots**2
                        if( sim_nuc_site_x(nuc_index) .ge. x_left(counter) .and. sim_nuc_site_x(nuc_index) .le. x_right(counter) .and. &
                                sim_nuc_site_z(nuc_index) .ge. z_bottom(counter) .and. sim_nuc_site_z(nuc_index) .le. z_top(counter)) then
                                waittime(nuc_index) = ht_tWait_hydrophobic
                       else
                        end if
                end do
            end do


!--------------------------------------------------------------------------------------------------
  do lb = 1,blockCount

         blockID = blockList(lb)
 
        call Grid_getDeltas(blockID,del)

        ! Get Blocks internal limits indexes:
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

        ! Point to blocks center and face vars:
        call Grid_getBlkPtr(blockID,solnData,CENTER)
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)
        call Grid_getBlkPtr(blockID,facezData,FACEZ)

        !---------------------------------------------------------------------
        !- mslee -  get interface velocity for level set advection

#if NDIM == 2

        call mph_getInterfaceVelocity( facexData(VELC_FACE_VAR,:,:,:),            &
                      faceyData(VELC_FACE_VAR,:,:,:),            &
                      ins_invRe,                                 &
                      blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                      blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                      del(DIR_X),del(DIR_Y), &
                      solnData(VISC_VAR,:,:,:), &
                      facexData(RH1F_FACE_VAR,:,:,:),            &
                      facexData(RH2F_FACE_VAR,:,:,:),            &
                      faceyData(RH1F_FACE_VAR,:,:,:),            &
                      faceyData(RH2F_FACE_VAR,:,:,:),            &
                      ins_gravX, ins_gravY, &
                      solnData(MDOT_VAR,:,:,:),&
                      solnData(SMRH_VAR,:,:,:),&
                      solnData(NRMX_VAR,:,:,:),&
                      solnData(NRMY_VAR,:,:,:),&
                      facexData(VELI_FACE_VAR,:,:,:),            &
                      faceyData(VELI_FACE_VAR,:,:,:),solnData(CURV_VAR,:,:,:))

#endif

#if NDIM == 3

        call mph_getInterfaceVelocity_3D( facexData(VELC_FACE_VAR,:,:,:),            &
                      faceyData(VELC_FACE_VAR,:,:,:),            &
                      facezData(VELC_FACE_VAR,:,:,:),            &
                      ins_invRe,                                 &
                      blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                      blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                      blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
                      del(DIR_X),del(DIR_Y),del(DIR_Z),&
                      solnData(VISC_VAR,:,:,:), &
                      facexData(RH1F_FACE_VAR,:,:,:),            &
                      facexData(RH2F_FACE_VAR,:,:,:),            &
                      faceyData(RH1F_FACE_VAR,:,:,:),            &
                      faceyData(RH2F_FACE_VAR,:,:,:),            &
                      facezData(RH1F_FACE_VAR,:,:,:),            &
                      facezData(RH2F_FACE_VAR,:,:,:),            &
                      ins_gravX, ins_gravY, ins_gravZ,&
                      solnData(MDOT_VAR,:,:,:),&
                      solnData(SMRH_VAR,:,:,:),&
                      solnData(NRMX_VAR,:,:,:),&
                      solnData(NRMY_VAR,:,:,:),&
                      solnData(NRMZ_VAR,:,:,:),&
                      facexData(VELI_FACE_VAR,:,:,:),            &
                      faceyData(VELI_FACE_VAR,:,:,:),            &
                      facezData(VELI_FACE_VAR,:,:,:),solnData(CURV_VAR,:,:,:))

#endif
         !---------------------------------------------------------------------


         ! Release pointers:
         call Grid_releaseBlkPtr(blockID,solnData,CENTER)
         call Grid_releaseBlkPtr(blockID,facexData,FACEX)
         call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
         call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  end do
 
  gcMask = .FALSE.
  gcMask(NUNK_VARS+VELI_FACE_VAR) = .TRUE.                 ! mslee - interface velocity
  gcMask(NUNK_VARS+1*NFACE_VARS+VELI_FACE_VAR) = .TRUE.    ! mslee - interface velocity

#if NDIM == 3
  gcMask(NUNK_VARS+2*NFACE_VARS+VELI_FACE_VAR) = .TRUE.
#endif

  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)

!---- Advancing/Receding Contact Angle -------!

do lb = 1,blockCount

     blockID = blockList(lb)

     call Grid_getBlkBoundBox(blockId,boundBox)

     bsize(:) = boundBox(2,:) - boundBox(1,:)

     call Grid_getBlkCenterCoords(blockId,coord)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     ! Point to blocks center and face vars:
        call Grid_getBlkPtr(blockID,solnData,CENTER)
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)

#if NDIM == 3
        call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif
!----------------------------------------------------------------------------------------Shantanu
        ! mph_psi(:,:,blockID)     = ht_psi
          mph_psi(:,:,blockID)     = sim_theta_hydrophilic
!-----------------------------------------------------------------------------------
         mph_psi_adv              = (90.0/180.0)*acos(-1.0)
         mph_vlim                 = 0.2
!-----------------------------------------------------------------------------------------------------------Shantanu
         do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
         do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

   !         if(solnData(DFUN_VAR,i,NGUARD+1,k)*solnData(DFUN_VAR,i+1,NGUARD+1,k) .le. 0 .or. &
   !            solnData(DFUN_VAR,i,NGUARD+1,k)*solnData(DFUN_VAR,i-1,NGUARD+1,k) .le. 0 .or. &
   !            solnData(DFUN_VAR,i,NGUARD+1,k)*solnData(DFUN_VAR,i,NGUARD+1,k+1) .le. 0 .or. &
   !            solnData(DFUN_VAR,i,NGUARD+1,k)*solnData(DFUN_VAR,i,NGUARD+1,k-1) .le. 0) then
   !
   !              veli = (facexData(VELI_FACE_VAR,i,NGUARD+1,k)+facexData(VELI_FACE_VAR,i+1,NGUARD+1,k))*0.5*&
   !                                solnData(NRMX_VAR,i,NGUARD+1,k) + &
   !                     (facezData(VELI_FACE_VAR,i,NGUARD+1,k)+facezData(VELI_FACE_VAR,i,NGUARD+1,k+1))*0.5*&
   !                                solnData(NRMZ_VAR,i,NGUARD+1,k)
   !             
   !              if(veli .ge. 0.0) then
   !              if(abs(veli) .le. mph_vlim) then
   !
   !                   mph_psi(i,k,blockID) = ((mph_psi_adv - ht_psi)/(2*mph_vlim))*abs(veli) + &
   !                                           (mph_psi_adv + ht_psi)/2.0d0
   !
   !              else
   !
   !                   mph_psi(i,k,blockID) = mph_psi_adv
   !
   !              end if
   !              end if
   !
   !         end if

                    x_position  = coord(IAXIS) - bsize(IAXIS)/2.0 +  &
                     		  real(i - NGUARD - 1)*del(IAXIS)  +  &
                    	          0.5*del(IAXIS)
        		    z_position  = coord(KAXIS) - bsize(KAXIS)/2.0 +  &
                                  real(k - NGUARD - 1)*del(KAXIS)  +  &
                                  0.5*del(KAXIS)
		    do nuc_index = 1,sim_nucSiteDens
                    	do counter = 1,nspots**2
                        	if( x_position .ge. x_left(counter) .and. x_position .le. x_right(counter) .and. &
                                    z_position .ge. z_bottom(counter) .and. z_position .le. z_top(counter)) then
                                   	 mph_psi(i,k,blockID) = sim_theta_hydrophobic
                        	else
                        	end if
                	end do
            	   end do


         end do
         end do
!---------------------------------------------------------------------------------------------------------------------------------
     ! Release pointers:
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM ==3
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

enddo

    !------------------------------------------------------------------
    !- kpd - Advect the multiphase distance function using WENO3 scheme
    !------------------------------------------------------------------

    call cpu_time(t_startMP2)

    volSum = 0.0
    volSumAll = 0.0

  do ii=1,1

    do lb = 1,blockCount

     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     ! Point to blocks center and face vars:
        call Grid_getBlkPtr(blockID,solnData,CENTER)
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)

#if NDIM == 3

        call Grid_getBlkPtr(blockID,facezData,FACEZ)

        !-----------------------------------------------------------------
        !Store Phi at previous time step for RK2
        if (ii.eq.1) solnData(AAJUNK_VAR,:,:,:) = solnData(DFUN_VAR,:,:,:)
        !-----------------------------------------------------------------

        !------------------------------------
        !! Call DFUN advection routine for 3D:
        !------------------------------------

       ! put ins_alfa*dt for RK2 ~~ Akash

        call mph_KPDadvectWENO3_3D(solnData(DFUN_VAR,:,:,:), &
                          facexData(VELI_FACE_VAR,:,:,:), &
                          faceyData(VELI_FACE_VAR,:,:,:), &
                          facezData(VELI_FACE_VAR,:,:,:), &
                          dt, &
                          del(DIR_X), &
                          del(DIR_Y), &
                          del(DIR_Z), &
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS), &
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS), &
                          blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),blockID)

        !KPD - Compute the Bubble Volume
        do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
              do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
                 if (solnData(DFUN_VAR,i,j,k) .gt. 0) then
                   volSum = volSum + (del(DIR_X) * del(DIR_Y) * del(DIR_Z))
                 end if
              end do
           end do
        end do

#elif NDIM == 2
        !-----------------------------------------------------------------
        !Store Phi at previous time step for RK2
        if (ii.eq.1) solnData(AAJUNK_VAR,:,:,:) = solnData(DFUN_VAR,:,:,:)
        !-----------------------------------------------------------------

        !------------------------------------
        ! Call DFUN advection routine for 2D:
        !------------------------------------

        !if ( lb .eq. 1 .AND. ins_meshMe .eq. 0) then
        !   print*,"ins_alfa, dt", ins_alfa, dt, del, blkLimits
        !endif

        ! put ins_alfa*dt for RK2 ~~ Akash

        call mph_KPDadvectWENO3(solnData(DFUN_VAR,:,:,:), &
                          facexData(VELI_FACE_VAR,:,:,:), &
                          faceyData(VELI_FACE_VAR,:,:,:), &
                          dt, &
                          del(DIR_X), &
                          del(DIR_Y), &
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),blockID)

        if(ii == 1) then
        !KPD - Compute the Bubble Volume
        do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
              if (solnData(DFUN_VAR,i,j,1) .gt. 0) then
                volSum = volSum + (del(DIR_X) * del(DIR_Y))
              end if
           end do
        end do
        endif

#endif

     ! Release pointers:
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM ==3
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
    enddo

    if(ii == 1) then
    call MPI_Allreduce(volSum, volSumAll, 1, FLASH_REAL,&
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    if (mph_meshMe .eq. 0) print*,"----------------------------------------"
    if (mph_meshMe .eq. 0) print*,"Total Liquid Volume: ",volSumAll
    if (mph_meshMe .eq. 0) print*,"----------------------------------------"
    endif

#if NDIM == 3

    mph_radius =  ((3.0*volSumAll)/(4*acos(-1.0)))**(1.0/3.0)

#endif

#if NDIM == 2

    mph_radius = sqrt(volSumAll/acos(-1.0))

#endif

    !********************************************************************************************************
    !-kpd - Fill distance function guard cells before re-initialization to
    !communicate updates
    gcMask = .FALSE.
    gcMask(DFUN_VAR) = .TRUE.
    gcMask(AAJUNK_VAR) = .TRUE.
#ifdef FLASH_GRID_PARAMESH
    intval = 1
    !intval = 2
    interp_mask_unk = intval;   interp_mask_unk_res = intval;
    interp_mask_work = intval;
#endif
    call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)
    !********************************************************************************************************
    !call mph_bcLevelSet(0.0d0,0.0d0)

    if (ii.eq.2) then
    do lb = 1,blockCount        ! Shizhao adds the loop over block
     blockID = blockList(lb)
       ! Point to blocks center and face vars:
       call Grid_getBlkPtr(blockID,solnData,CENTER)

       solnData(DFUN_VAR,:,:,:) = 0.5* (solnData(AAJUNK_VAR,:,:,:) + solnData(DFUN_VAR,:,:,:))

       ! Release pointers:
       call Grid_releaseBlkPtr(blockID,solnData,CENTER)
    enddo
    end if
 
 end do   !End ii RK loop


    !********************************************************************************************************
    !-kpd - Fill distance function guard cells before re-initialization to
    !communicate updates
    gcMask = .FALSE.
    gcMask(DFUN_VAR) = .TRUE.
#ifdef FLASH_GRID_PARAMESH
    intval = 1
    !intval = 2
    interp_mask_unk = intval;   interp_mask_unk_res = intval;
    interp_mask_work = intval;
#endif
    call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)
    !********************************************************************************************************
    !call mph_bcLevelSet(0.0d0,0.0d0)

    !call sim_trackBubble(blockCount, blockList,vol, cx, cy, vx, vy,xh,yh,xl,yl)
    !if(mph_meshMe == 0) write(100000,'(12f20.12)') real(dr_nStep), timeEndAdv,vol, cx, cy, vx, vy, xh, yh, xl, yl

!###################################################################################################################

! Level set re-initialization for nucleate boiling simulation.
! The procedure checks if the bubble has departed the heating surface (boundary),
! and if it has then it creates a nucleus for the next bubble to grow and
! continue the cyclic process.
! - Akash Dhruv

#ifdef NUCLEATE_BOILING

if(ins_meshMe .eq. MASTER_PE) then
        print *,"Nucleation site, truth value, temperature and time stamp"
        do nuc_index = 1,sim_nucSiteDens
              print *,mph_isAttachedAll(nuc_index),mph_nucSiteTemp(nuc_index),mph_timeStampAll(nuc_index)
        end do
end if

!if(ins_meshMe .eq. MASTER_PE)write (*,*)"Nucleation site truth value - ",mph_isAttachedAll(1:sim_nucSiteDens)
!if(ins_meshMe .eq. MASTER_PE)write (*,*)"Nucleation site temperature - ",(mph_nucSiteTemp(1:sim_nucSiteDens))
!if(ins_meshMe .eq. MASTER_PE)write (*,*)"Nucleation site time        - ",(mph_timeStampAll(1:sim_nucSiteDens))

do nuc_index =1,sim_nucSiteDens

  isAttached  = .false.
  nucSiteTemp = 0.0

  do lb = 1,blockCount

     blockID = blockList(lb)
     call Grid_getBlkBoundBox(blockId,boundBox)

     bsize(:) = boundBox(2,:) - boundBox(1,:)

     call Grid_getBlkCenterCoords(blockId,coord)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)

     ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
              real(blkLimits(LOW,JAXIS) - NGUARD - 1)*del(JAXIS)  +  &
              0.5*del(JAXIS)

     if(abs(ycell-0.5*del(JAXIS)) .le. tol) then

       do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)-1
        do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)-1

           xcell  = coord(IAXIS) - bsize(IAXIS)/2.0 +  &
                    real(i - NGUARD - 1)*del(IAXIS)  +  &
                    0.5*del(IAXIS)

           xcellp = coord(IAXIS) - bsize(IAXIS)/2.0 +  &
                    real(i+1 - NGUARD - 1)*del(IAXIS)  +  &
                    0.5*del(IAXIS)

           zcell  = coord(KAXIS) - bsize(KAXIS)/2.0 +  &
                    real(k - NGUARD - 1)*del(KAXIS)  +  &
                    0.5*del(KAXIS)

           zcellp = coord(KAXIS) - bsize(KAXIS)/2.0 +  &
                    real(k+1 - NGUARD - 1)*del(KAXIS)  +  &
                    0.5*del(KAXIS)

           if((xcell .le. sim_nuc_site_x(nuc_index)) .and. (xcellp .ge. sim_nuc_site_x(nuc_index)) .and. &
              (zcell .le. sim_nuc_site_z(nuc_index)) .and. (zcellp .ge. sim_nuc_site_z(nuc_index))) then

             if ((solnData(DFUN_VAR,i,blkLimits(LOW,JAXIS),k)   .ge. 0.0) .or. (solnData(DFUN_VAR,i+1,blkLimits(LOW,JAXIS),k)   .ge. 0.0) .or. &
                 (solnData(DFUN_VAR,i,blkLimits(LOW,JAXIS),k+1) .ge. 0.0) .or. (solnData(DFUN_VAR,i+1,blkLimits(LOW,JAXIS),k+1) .ge. 0.0)) then
                        
                        isAttached = isAttached .or. .true.
             else
                        isAttached = isAttached .or. .false.
             end if 

             nucSiteTemp = (solnData(TEMP_VAR,i,blkLimits(LOW,JAXIS),k) + &
                            solnData(TEMP_VAR,i+1,blkLimits(LOW,JAXIS),k) + &
                            solnData(TEMP_VAR,i,blkLimits(LOW,JAXIS),k+1) + & 
                            solnData(TEMP_VAR,i+1,blkLimits(LOW,JAXIS),k+1))/4.0

           end if

        end do
       end do

     end if
  
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  end do

  call MPI_Allreduce(isAttached, mph_isAttachedAll(nuc_index), 1, FLASH_LOGICAL,&
                     MPI_LOR, MPI_COMM_WORLD, ierr)

  call MPI_Allreduce(nucSiteTemp, mph_nucSiteTemp(nuc_index), 1, FLASH_REAL,&
                     MPI_MAX, MPI_COMM_WORLD, ierr)

  if((mph_isAttachedOld(nuc_index) .eqv. .true.) .and. (mph_isAttachedAll(nuc_index) .eqv. .false.)) then

         mph_timeStampAll(nuc_index) = dr_simTime
         !mph_offset = mph_offset+mph_radius

  end if

  mph_isAttachedOld(nuc_index) = mph_isAttachedAll(nuc_index)
!-----------------------------------------------------------------------------------------Shantanu
  if( (mph_isAttachedAll(nuc_index) .eqv. .false.) .and. &
      (mph_timeStampAll(nuc_index) + waittime(nuc_index) .le. dr_simTime) .and. &
      (mph_nucSiteTemp(nuc_index) .ge. ht_Tnuc) )then
!  if( (mph_isAttachedAll(nuc_index) .eqv. .false.) .and. &
!      (mph_timeStampAll(nuc_index) +  0.2 .le. dr_simTime) .and. &
!      (mph_nucSiteTemp(nuc_index) .ge. ht_Tnuc) )then

!-----------------------------------------------------------------------------------------------------
  do lb = 1,blockCount

     blockID = blockList(lb)
     call Grid_getBlkBoundBox(blockId,boundBox)

     bsize(:) = boundBox(2,:) - boundBox(1,:)

     call Grid_getBlkCenterCoords(blockId,coord)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)

    do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
      do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

         xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
                 real(i - NGUARD - 1)*del(IAXIS) +   &
                 0.5*del(IAXIS)

         ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                  real(j - NGUARD - 1)*del(JAXIS)  +  &
                  0.5*del(JAXIS)

         zcell  = coord(KAXIS) - bsize(KAXIS)/2.0 + &
                  real(k - NGUARD - 1)*del(KAXIS)  + &
                  0.5*del(KAXIS)

         nuc_dfun  = 0.1 - sqrt((xcell-sim_nuc_site_x(nuc_index))**2+(ycell-sim_nuc_site_y(nuc_index))**2+(zcell-sim_nuc_site_z(nuc_index))**2)

         solnData(DFUN_VAR,i,j,k) = max(solnData(DFUN_VAR,i,j,k),nuc_dfun)

      end do
     end do
    end do

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  end do

    gcMask = .FALSE.
    gcMask(DFUN_VAR) = .TRUE.
#ifdef FLASH_GRID_PARAMESH
    intval = 1
    !intval = 2
    interp_mask_unk = intval;   interp_mask_unk_res = intval;
    interp_mask_work = intval;
#endif
    call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)

  end if

end do

#endif
! End of procedure - Akash
!###################################################################################################################

   call cpu_time(t_startMP2a)

   do ii = 1,mph_lsit

     !------------------------------
     !- kpd - Level set redistancing 
     !------------------------------

     !lsDT = dt
     lsT  = 0.0

     do lb = 1,blockCount
        blockID = blockList(lb)

         call Grid_getBlkBoundBox(blockId,boundBox)
         bsize(:) = boundBox(2,:) - boundBox(1,:)

        call Grid_getBlkCenterCoords(blockId,coord)

        ! Get blocks dx, dy ,dz:
        call Grid_getDeltas(blockID,del)

        ! Get Blocks internal limits indexes:
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

        ! Point to blocks center and face vars:
        call Grid_getBlkPtr(blockID,solnData,CENTER)
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
        call Grid_getBlkPtr(blockID,facezData,FACEZ)

        !--------------------------------------------
        ! Call DFUN re-initialization routine for 3D:
        !--------------------------------------------
        lsDT = MIN(10.0*dt,0.001)
        !minCellDiag = SQRT(del(DIR_X)**2.+del(DIR_Y)**2.+del(DIR_Z)**2.)
        minCellDiag =SQRT((SQRT(del(DIR_X)**2.+del(DIR_Y)**2.))**2.+del(DIR_Z)**2.)
        !lsDT = minCellDiag/2.0d0
        if ( ii .eq. mph_lsit .AND. lb .eq. 1 .AND. mph_meshMe .eq. 0) then
           print*,"Level Set Initialization Iteration # ",ii,minCellDiag,lsDT
        end if

        if (ii.eq.1) solnData(AAJUNK_VAR,:,:,:) = solnData(DFUN_VAR,:,:,:)

        call mph_KPDlsRedistance_3D(solnData(DFUN_VAR,:,:,:), &
                          facexData(VELI_FACE_VAR,:,:,:), &
                          faceyData(VELI_FACE_VAR,:,:,:), &
                          facezData(VELI_FACE_VAR,:,:,:), &
                          del(DIR_X),del(DIR_Y),del(DIR_Z),  &
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS), &
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS), &
                          blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS), &
                          solnData(AAJUNK_VAR,:,:,:), lsDT, minCellDiag )

#elif NDIM == 2

        !--------------------------------------------
        ! Call DFUN re-initialization routine for 2D:
        !--------------------------------------------
        !lsDT = MIN(10.0*dt,0.001)
        minCellDiag = SQRT(del(DIR_X)**2.+del(DIR_Y)**2.)
        lsDT = minCellDiag/2.0d0
        if ( ii .eq. mph_lsit .AND. lb .eq. 1 .AND. mph_meshMe .eq. 0) then
           print*,"Level Set Initialization Iteration # ",ii,minCellDiag,lsDT
        end if

        if (ii.eq.1) solnData(AAJUNK_VAR,:,:,:) = solnData(DFUN_VAR,:,:,:)

        call mph_KPDlsRedistance(solnData(DFUN_VAR,:,:,:), &
                          facexData(VELI_FACE_VAR,:,:,:),  &
                          faceyData(VELI_FACE_VAR,:,:,:),  &
                          del(DIR_X),del(DIR_Y),  &
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS), &
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS), &
                          solnData(AAJUNK_VAR,:,:,:), lsDT, blockID,minCellDiag)

          !if(ii .eq. mph_lsit) then

         !print *,"DFUN - ",solnData(DFUN_VAR,NXB/2,NYB/2,1)," - PROC - ",mph_meshME

        !end if

#endif

        ! Release pointers:
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM ==3
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif
     enddo

    !*********************************************************************************************************
    !-kpd - Fill distance function guard cells after each re-initialization to
    !communicate updates
    gcMask = .FALSE.
    gcMask(DFUN_VAR) = .TRUE.
#ifdef FLASH_GRID_PARAMESH
    intval = 1
    !intval = 2
    interp_mask_unk = intval;   interp_mask_unk_res = intval;
#endif
    call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)
    !*********************************************************************************************************
    !call mph_bcLevelSet(0.0d0,0.0d0)

!***********************************************************************************************
!***********************************************************************************************
!***********************************************************************************************

      lsT = lsT + lsDT

   end do  ! End do: ii=1,lsit

   call cpu_time(t_stopMP2)
   if (mph_meshMe .eq. 0) print*,"Total Multiphase Time: ",t_stopMP2-t_startMP2,t_stopMP2-t_startMP2a

   !print*,"Multiphase 2 Solver Time  ",t_stopMP2-t_startMP2

end subroutine
