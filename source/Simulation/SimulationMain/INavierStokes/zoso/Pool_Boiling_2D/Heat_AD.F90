subroutine Heat_AD( blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)

! This subroutine solves the heat advection diffusion equation, calculates
! heat and mass flux for phase change simulations

! Author - A.V. Dhruv

#include "constants.h"
#include "Heat_AD.h"
#include "Flash.h"

   use Heat_AD_interface, only: Heat_Solve,Heat_RHS_upwind,Heat_calGradT,Heat_calGradT_central,&
                                Heat_extrapGradT,Heat_calMdot,Heat_RHS_3D,Heat_RHS_weno3,&
                                Heat_extrapGradT_3D,Heat_calGradT_3D,Heat_RHS_central,&
                                Heat_RHS_3D_weno3,Heat_calGradT_3D_central,Heat_extrapGradT_weno3,Heat_getQmicro,&
                                Heat_getWallflux

   use Grid_interface, only: Grid_getDeltas, Grid_getBlkIndexLimits,&
                             Grid_getBlkPtr, Grid_releaseBlkPtr,    &
                             Grid_fillGuardCells 

   use IncompNS_data, only:ins_invRe,ins_meshMe

   use Multiphase_data, only: mph_thco1,mph_cp1,mph_thco2,mph_cp2,mph_meshMe,mph_meshNumProcs,mph_rho2,mph_baseRadius,mph_baseCountAll,mph_baseCount

   use Heat_AD_data, only: ht_hfit,ht_Tsat,ht_microFlg,ht_fmic,ht_qmic

   use Driver_data,  only: dr_nstep,dr_simTime

   use Heat_AD_data, only: ht_AMR_specs, ht_qmic, ht_dxmin, ht_Nu_l, ht_Nu_t

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

   implicit none

   include "Flash_mpi.h"

   integer, INTENT(INOUT) :: blockCount
   integer, INTENT(IN) :: sweepOrder
   integer, INTENT(INOUT) :: blockList(MAXBLOCKS)
   real,    INTENT(IN) :: timeEndAdv, dt, dtOld

   integer ::  blockID,lb,step,ierr,i,j,k
   real :: xcell,ycell
   real ::  del(MDIM),coord(MDIM),bsize(MDIM)
   real ::  boundBox(2,MDIM)
   integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
   logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)
   real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData
   real :: Tnl_res, Tnv_res, Tnl_res1, Tnv_res1,Tnl_resBlock,Tnv_resBlock
   real :: T_res,T_res1,T_resBlock
   real :: max_flux, min_flux
   real :: Nu_l,Nu_t
   integer :: hcounterAll,hcounter
   integer :: iter_count,intval
   real    :: dxmin

   iter_count = 0

!_________________________________Microlayer solution_____________________________________________!

   if(ht_microFlg) then

      print *,"Calculating minimum dx: "

      dxmin    = 1e10

      do lb = 1,blockCount

        blockID = blockList(lb)
        call Grid_getDeltas(blockID,del)
        dxmin = min(dxmin,del(JAXIS))

      end do

      call MPI_ALLREDUCE(dxmin,ht_dxmin,1,FLASH_REAL,MPI_MIN,MPI_COMM_WORLD,ierr)
       
      ht_microFlg = .FALSE. 

   end if

   !if (ins_meshMe .eq. MASTER_PE) call Heat_getQmicro(ht_qmic,ht_fmic,ht_dxmin)

   !call MPI_BCAST(ht_qmic, 1, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
   !call MPI_BCAST(ht_fmic, 1, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

   if (ins_meshMe .eq. MASTER_PE) print *,"qmic,fmic: ",ht_qmic,ht_fmic

!______________________________________Energy Equation____________________________________________!

   T_resBlock   = 0.0

   !if (dr_simTime .ge. 1600.00 .and. dr_simTime .le. 1900.00) then

   !  ht_Tsat  = 0.0013*(dr_simTime-1600.00) + 0.0
   !  mph_rho2 = 170 - 0.1*(dr_simTime-1600.00)

   !  if (ins_meshMe .eq. MASTER_PE) call Heat_getQmicro(ht_qmic,ht_dxmin)

   !  call MPI_BCAST(ht_qmic, 1, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

   !  print *,"qmic: ",ht_qmic

   !end if

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
     call Grid_getBlkPtr(blockID,facezData,FACEZ)

     ! Calculate RHS for advections diffusion
#if NDIM == 2
     call Heat_RHS_weno3(solnData(RHST_VAR,:,:,:), solnData(TEMP_VAR,:,:,:),&
                     facexData(VELC_FACE_VAR,:,:,:),&
                     faceyData(VELC_FACE_VAR,:,:,:),&
                     del(DIR_X),del(DIR_Y),del(DIR_Z),&
                     ins_invRe,&
                     blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                     blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                     facexData(RH1F_FACE_VAR,:,:,:),facexData(RH2F_FACE_VAR,:,:,:),&
                     faceyData(RH1F_FACE_VAR,:,:,:),faceyData(RH2F_FACE_VAR,:,:,:),&
                     solnData(ALPH_VAR,:,:,:),&
                     solnData(PFUN_VAR,:,:,:),solnData(DFUN_VAR,:,:,:),&
                     solnData(MDOT_VAR,:,:,:),solnData(NRMX_VAR,:,:,:),&
                     solnData(NRMY_VAR,:,:,:),solnData(SMRH_VAR,:,:,:),&
                     solnData(CURV_VAR,:,:,:))
#endif

#if NDIM == 3
     call Heat_RHS_3D_weno3(solnData(RHST_VAR,:,:,:), solnData(TEMP_VAR,:,:,:),&
                     facexData(VELC_FACE_VAR,:,:,:),&
                     faceyData(VELC_FACE_VAR,:,:,:),&
                     facezData(VELC_FACE_VAR,:,:,:),&
                     del(DIR_X),del(DIR_Y),del(DIR_Z),&
                     ins_invRe,&
                     blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                     blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                     blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
                     facexData(RH1F_FACE_VAR,:,:,:),facexData(RH2F_FACE_VAR,:,:,:),&
                     faceyData(RH1F_FACE_VAR,:,:,:),faceyData(RH2F_FACE_VAR,:,:,:),&
                     facezData(RH1F_FACE_VAR,:,:,:),facezData(RH2F_FACE_VAR,:,:,:),&
                     solnData(ALPH_VAR,:,:,:),&
                     solnData(PFUN_VAR,:,:,:),solnData(DFUN_VAR,:,:,:),&
                     solnData(MDOT_VAR,:,:,:),solnData(NRMX_VAR,:,:,:),&
                     solnData(NRMY_VAR,:,:,:),solnData(NRMZ_VAR,:,:,:),&
                     solnData(SMRH_VAR,:,:,:),solnData(CURV_VAR,:,:,:))
#endif

     ! RK-1 dt/1.0
     ! RK-2 dt/2.0

     ! Calculate temperature at new time-step

     call Heat_Solve(solnData(TEMP_VAR,:,:,:),solnData(RHSO_VAR,:,:,:),&
                     solnData(RHST_VAR,:,:,:),&
                     dt,&
                     blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                     blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                     blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),T_res1)


     ! uncomment for RK-1
     T_resBlock = T_resBlock + T_res1

     solnData(TNLQ_VAR,:,:,:) = 0.0
     solnData(TNVP_VAR,:,:,:) = 0.0

     solnData(RHSO_VAR,:,:,:) = solnData(RHST_VAR,:,:,:)

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

    end do

    gcMask = .FALSE.
    gcMask(TEMP_VAR)=.TRUE.

    call Grid_fillGuardCells(CENTER,ALLDIR,&
         maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask,selectBlockType=ACTIVE_BLKS)

   if(blockCount .gt. 0) T_resBlock = T_resBlock/blockCount

   ! Collect residuals from other processes
   call MPI_Allreduce(T_resBlock, T_res, 1, FLASH_REAL,&
                      MPI_SUM, MPI_COMM_WORLD, ierr)

   T_res = sqrt(T_res/mph_meshNumProcs)

   if(mph_meshMe .eq. MASTER_PE) print *,"T_res:",T_res

!___________________________________End of Energy Equation____________________________________________!

!____________________________________Heat Flux calculation____________________________________________!

   if(mph_meshMe .eq. MASTER_PE) print *,"Entering heat flux calculation"

   Nu_l = 0.0
   Nu_t = 0.0
   hcounter = 0

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
    
     ! Calculate Heat Flux in known phase

#if NDIM == 2
     call Heat_calGradT_central(solnData(TNLQ_VAR,:,:,:),solnData(TNVP_VAR,:,:,:),&
                        solnData(TEMP_VAR,:,:,:),solnData(DFUN_VAR,:,:,:),&
                        solnData(PFUN_VAR,:,:,:),del(DIR_X),del(DIR_Y),del(DIR_Z),&
                        blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                        blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                        solnData(NRMX_VAR,:,:,:),solnData(NRMY_VAR,:,:,:),solnData(MFLG_VAR,:,:,:))
#endif

#if NDIM == 3
     call Heat_calGradT_3D_central(solnData(TNLQ_VAR,:,:,:),solnData(TNVP_VAR,:,:,:),&
                        solnData(TEMP_VAR,:,:,:),solnData(DFUN_VAR,:,:,:),&
                        solnData(PFUN_VAR,:,:,:),del(DIR_X),del(DIR_Y),del(DIR_Z),&
                        blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                        blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                        blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
                        solnData(NRMX_VAR,:,:,:),solnData(NRMY_VAR,:,:,:),solnData(NRMZ_VAR,:,:,:),&
                        solnData(MFLG_VAR,:,:,:))
#endif

     ! Microlayer contribution - only for nucleate boiling
     !if(dr_nstep .gt. 1) then

     solnData(TNLQ_VAR,:,:,:) = solnData(TNLQ_VAR,:,:,:) + solnData(TMIC_VAR,:,:,:)*(ht_qmic/(2.0*del(DIR_X)*del(DIR_Y)))

     !solnData(TNLQ_VAR,:,:,:) = solnData(TNLQ_VAR,:,:,:) + solnData(TMIC_VAR,:,:,:)*(ht_qmic/(del(DIR_X)))

     !   solnData(RTES_VAR,:,:,:) = solnData(TMIC_VAR,:,:,:)*((ht_qmic*mph_baseRadius)/(mph_baseCountAll*del(DIR_X)*del(DIR_Z)))

     !end if

     ! Release pointers
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  
   end do

   ! Apply BC
   gcMask = .FALSE.
   gcMask(TNLQ_VAR)=.TRUE.
   gcMask(TNVP_VAR)=.TRUE.

   call Grid_fillGuardCells(CENTER,ALLDIR,&
        maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask,selectBlockType=ACTIVE_BLKS)

   do lb=1,blockCount

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
 
     ! Wall heat flux

     ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
              real(blkLimits(LOW,JAXIS) - NGUARD - 1)*del(JAXIS)  +  &
              0.5*del(JAXIS)

     call Heat_getWallflux(solnData(PFUN_VAR,:,:,:),solnData(TEMP_VAR,:,:,:),solnData(RTES_VAR,:,:,:),Nu_l,Nu_t,hcounter,del(JAXIS),ycell,&
                          blkLimits(LOW,JAXIS),blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),blockID)

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
 
   end do

   call MPI_Allreduce(hcounter, hcounterAll, 1, FLASH_INTEGER,&
                      MPI_SUM, MPI_COMM_WORLD, ierr)

   call MPI_Allreduce(Nu_l, ht_Nu_l, 1, FLASH_REAL,&
                      MPI_SUM, MPI_COMM_WORLD, ierr)

   call MPI_Allreduce(Nu_t, ht_Nu_t, 1, FLASH_REAL,&
                      MPI_SUM, MPI_COMM_WORLD, ierr)

   ht_Nu_l = ht_Nu_l/hcounterAll
   ht_Nu_t = ht_Nu_t/hcounterAll

   if(ins_meshMe .eq. MASTER_PE) print *,"Wall Nusselt Number Liq - ",ht_Nu_l
   if(ins_meshMe .eq. MASTER_PE) Print *,"Wall Nusselt Number Tot - ",ht_Nu_t
  
   ! Apply BC
   gcMask = .FALSE.
   gcMask(RTES_VAR)=.TRUE.

   call Grid_fillGuardCells(CENTER,ALLDIR,&
        maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask,selectBlockType=ACTIVE_BLKS) 
!_________________________________End of Heat Flux calculation_____________________________!

!__________________________Heat Flux extrapolation sub iterations__________________________!

   do while(iter_count <= ht_hfit)

   Tnl_resBlock = 0.0
   Tnv_resBlock = 0.0

   do lb = 1,blockCount

     Tnl_res1 = 0.0
     Tnv_res1 = 0.0

     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     ! Lock pointers
     call Grid_getBlkPtr(blockID,solnData,CENTER)

     ! Extrapolation of heat flux into unknown phase

#if NDIM == 2
     call Heat_extrapGradT(solnData(TNLQ_VAR,:,:,:),solnData(TNVP_VAR,:,:,:),&
                           solnData(TEMP_VAR,:,:,:),solnData(DFUN_VAR,:,:,:),&
                           solnData(PFUN_VAR,:,:,:),del(DIR_X),del(DIR_Y),del(DIR_Z),&
                           solnData(NRMX_VAR,:,:,:),solnData(NRMY_VAR,:,:,:),&
                           blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                           blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                           Tnl_res1,Tnv_res1,solnData(MFLG_VAR,:,:,:))
#endif

#if NDIM == 3
     call Heat_extrapGradT_3D(solnData(TNLQ_VAR,:,:,:),solnData(TNVP_VAR,:,:,:),&
                           solnData(TEMP_VAR,:,:,:),solnData(DFUN_VAR,:,:,:),&
                           solnData(PFUN_VAR,:,:,:),del(DIR_X),del(DIR_Y),del(DIR_Z),&
                           solnData(NRMX_VAR,:,:,:),solnData(NRMY_VAR,:,:,:),&
                           solnData(NRMZ_VAR,:,:,:),&
                           blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                           blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                           blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
                           Tnl_res1,Tnv_res1,solnData(MFLG_VAR,:,:,:))
#endif
     ! Release pointers

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)

     Tnl_resBlock = Tnl_resBlock + Tnl_res1
     Tnv_resBlock = Tnv_resBlock + Tnv_res1
   
   end do

     if(blockCount .gt. 0) Tnl_resBlock = Tnl_resBlock/blockCount
     if(blockCount .gt. 0) Tnv_resBlock = Tnv_resBlock/blockCount

     ! Collect residuals from other processes
     call MPI_Allreduce(Tnl_resBlock, Tnl_res, 1, FLASH_REAL,&
                     MPI_SUM, MPI_COMM_WORLD, ierr)

     call MPI_Allreduce(Tnv_resBlock, Tnv_res, 1, FLASH_REAL,&
                     MPI_SUM, MPI_COMM_WORLD, ierr)

     Tnl_res = sqrt(Tnl_res/mph_meshNumProcs)
     Tnv_res = sqrt(Tnv_res/mph_meshNumProcs)

     ! Apply BC
     gcMask = .FALSE.
     gcMask(TNLQ_VAR)=.TRUE.
     gcMask(TNVP_VAR)=.TRUE.

     call Grid_fillGuardCells(CENTER,ALLDIR,&
         maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask,selectBlockType=ACTIVE_BLKS)

     ! Convergence check
     !if ((Tnl_res+Tnv_res)/2.d0 < 1E-6 ) exit
     if  (Tnl_res < 1E-6) exit

     ! Increment counter
     iter_count = iter_count + 1

   end do

     ! Print results
   if(mph_meshMe .eq. MASTER_PE) then

      print *,"Tnl_res:",Tnl_res
      print *,"Tnv_res:",Tnv_res
      print *,"Extrapolation Iterations (Max Iterations): ",iter_count,"(",ht_hfit,")"
     
   end if

!______End Heat Flux extrapolation sub iterations____________________________________!

!______Mass Flux calculation_________________________________________________________!

   max_flux = -10.**(10.)
   min_flux =  10.**(10.)

   do lb = 1,blockCount

     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     ! Lock pointers 
     call Grid_getBlkPtr(blockID,solnData,CENTER)

     ! Calculate Mass Flux
     call Heat_calMdot(solnData(MDOT_VAR,:,:,:),&
                       solnData(TNLQ_VAR,:,:,:),solnData(TNVP_VAR,:,:,:),&
                       mph_thco2,mph_thco1,&
                       solnData(NRMX_VAR,:,:,:),solnData(NRMY_VAR,:,:,:),&
                       blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                       blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                       blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS))

    ! Release pointers

    max_flux = max(max_flux,maxval(solnData(MDOT_VAR,:,:,:)))
    min_flux = min(min_flux,minval(solnData(MDOT_VAR,:,:,:)))

    call Grid_releaseBlkPtr(blockID,solnData,CENTER)
   
   end do

   ! Apply BC
   gcMask = .FALSE.
   gcMask(MDOT_VAR)=.TRUE.

   call Grid_fillGuardCells(CENTER,ALLDIR,&
        maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask,selectBlockType=ACTIVE_BLKS)

   call MPI_Allreduce(ht_AMR_specs(2), max_flux, 1, FLASH_REAL,&
                      MPI_MAX, MPI_COMM_WORLD, ierr)

   call MPI_Allreduce(ht_AMR_specs(1), min_flux, 1, FLASH_REAL,&
                      MPI_MIN, MPI_COMM_WORLD, ierr)

!______End Mass Flux calculation___________________________________________________!

end subroutine Heat_AD
