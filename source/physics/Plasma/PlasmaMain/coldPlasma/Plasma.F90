subroutine Plasma( blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)

   use Plasma_data
   use Plasma_interface, only: Plasma_Solve, Plasma_hvDiffCoeff,     &
                               Plasma_elDiffCoeff,Plasma_ColFreq,    &
                               Plasma_sumNeutrals,Plasma_spReactions,&
                               Plasma_spGeneration,Plasma_sumIons

   use Grid_interface, only: Grid_getDeltas, Grid_getBlkIndexLimits,&
                             Grid_getBlkPtr, Grid_releaseBlkPtr,    &
                             Grid_fillGuardCells
   implicit none
#include "constants.h"
#include "Plasma.h"
#include "Flash.h"   

!#define DEBUG_PLASMA

   include "Flash_mpi.h"

   integer, INTENT(INOUT) :: blockCount
   integer, INTENT(IN) :: sweepOrder
   integer, INTENT(INOUT) :: blockList(MAXBLOCKS)
   real,    INTENT(IN) :: timeEndAdv, dt, dtOld

   integer ::  blockID,lb
   real ::  del(MDIM)
   integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
   logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)
   real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData
   real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: oldT

   real :: T_res1,T_res,T_resBlock
   real, dimension(10) :: T_resHV, T_resBlockHV
   integer :: ierr, i,j,k

   T_resBlock      = 0.0
   T_resBlockHV(:) = 0.0

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
     
     !add all neutral species   
#ifdef DEBUG_PLASMA
     print *,"Going into sumNeurtral"
#endif
     solnData(DNAT_VAR,:,:,:) = 0.0
     do i = 0,5
       call Plasma_sumNeutrals( solnData(DHV0_VAR+i,:,:,:),solnData(DNAT_VAR,:,:,:), &
                                blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                                blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS))
     end do 
 
     !add all ions (+ and -)
#ifdef DEBUG_PLASMA
     print *,"Going into sumIons"
#endif
     solnData(DNIT_VAR,:,:,:) = 0.0
     do i = 0,3
       call Plasma_sumIons(solnData(DHV6_VAR+i,:,:,:),solnData(DNIT_VAR,:,:,:),&
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS))
     end do
 
     !obtain diffusion coefficient of neutral species
#ifdef DEBUG_PLASMA
     print *,"Going into hvDiffCoeff"
#endif
     do i = 0,5
       call Plasma_hvDiffCoeff(solnData(DFH0_VAR+i,:,:,:), &
                               solnData(PRHV_VAR,:,:,:),solnData(TPHV_VAR,:,:,:), &
                               blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                               blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                               pls_RSCD(i+1), pls_MHSP(i+1)) 
     end do
     
     !obtain collision frequencies
#ifdef DEBUG_PLASMA
     print *,"Going into ColFreq"
#endif
     call Plasma_ColFreq(solnData(FVEI_VAR,:,:,:), solnData(FVEA_VAR,:,:,:), &
                         solnData(DNAT_VAR,:,:,:), solnData(DELE_VAR,:,:,:), &
                         solnData(TPEL_VAR,:,:,:), &
                         blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                         blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS) )


     !obtain diffusion coefficient of electrons (and ions)
#ifdef DEBUG_PLASMA
     print *,"Going into elDiffCoeff"
#endif
     call Plasma_elDiffCoeff(solnData(DFIO_VAR,:,:,:),solnData(DFEL_VAR,:,:,:),&
                             solnData(TPEL_VAR,:,:,:),solnData(TPHV_VAR,:,:,:),&
                             solnData(FVEA_VAR,:,:,:),solnData(FVEI_VAR,:,:,:),&
                             blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                             blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS) )
      
     !obtain chemical reaction rates for all species
#ifdef DEBUG_PLASMA
     print *,"Going into spReactions"
#endif
     call Plasma_spReactions(solnData(RSP0_VAR,:,:,:), solnData(RSP1_VAR,:,:,:),&
                             solnData(RSP2_VAR,:,:,:), solnData(RSP3_VAR,:,:,:),&
                             solnData(RSP4_VAR,:,:,:), solnData(RSP5_VAR,:,:,:),&
                             solnData(RSP6_VAR,:,:,:), solnData(RSP7_VAR,:,:,:),&
                             solnData(RSP8_VAR,:,:,:), solnData(RSP9_VAR,:,:,:),&
                             solnData(RSP10_VAR,:,:,:), solnData(RSP11_VAR,:,:,:),&
                             solnData(RSP12_VAR,:,:,:), solnData(RSP13_VAR,:,:,:),&
                             solnData(RSP14_VAR,:,:,:), solnData(RSP15_VAR,:,:,:),&
                             solnData(RSP16_VAR,:,:,:),&
                             solnData(TPHV_VAR,:,:,:),  solnData(TPEL_VAR,:,:,:),&
                             blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                             blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS) )

     !obtain species generation rates using chemical reations
#ifdef DEBUG_PLASMA
     print *,"Going into spGeneration"
#endif
     call Plasma_spGeneration(solnData(DHV0_VAR,:,:,:), solnData(DHV1_VAR,:,:,:),& 
                              solnData(DHV2_VAR,:,:,:), solnData(DHV3_VAR,:,:,:),&
                              solnData(DHV4_VAR,:,:,:), solnData(DHV5_VAR,:,:,:),&
                              solnData(DHV6_VAR,:,:,:), solnData(DHV7_VAR,:,:,:),&
                              solnData(DHV8_VAR,:,:,:), solnData(DHV9_VAR,:,:,:),&
                              solnData(DELE_VAR,:,:,:),                          &
                              solnData(RSP0_VAR,:,:,:), solnData(RSP1_VAR,:,:,:),&
                              solnData(RSP2_VAR,:,:,:), solnData(RSP3_VAR,:,:,:),&
                              solnData(RSP4_VAR,:,:,:), solnData(RSP5_VAR,:,:,:),&
                              solnData(RSP6_VAR,:,:,:), solnData(RSP7_VAR,:,:,:),&
                              solnData(RSP8_VAR,:,:,:), solnData(RSP9_VAR,:,:,:),&
                              solnData(RSP10_VAR,:,:,:),solnData(RSP11_VAR,:,:,:),&
                              solnData(RSP12_VAR,:,:,:),solnData(RSP13_VAR,:,:,:),&
                              solnData(RSP14_VAR,:,:,:), solnData(RSP15_VAR,:,:,:),&
                              solnData(RSP16_VAR,:,:,:),&
                              solnData(GNH0_VAR,:,:,:), solnData(GNH1_VAR,:,:,:),& 
                              solnData(GNH2_VAR,:,:,:), solnData(GNH3_VAR,:,:,:),&
                              solnData(GNH4_VAR,:,:,:), solnData(GNH5_VAR,:,:,:),&
                              solnData(GNH6_VAR,:,:,:), solnData(GNH7_VAR,:,:,:),&
                              solnData(GNH8_VAR,:,:,:), solnData(GNH9_VAR,:,:,:),&
                              solnData(GNE_VAR,:,:,:),  solnData(GNEBZ_VAR,:,:,:),&
                              solnData(GNERT_VAR,:,:,:),                          &
                              blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                              blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS) )     
 
     !obtain number density of electrons
#ifdef DEBUG_PLASMA
     print *,"Going into Plasma solve electrons"
#endif
     oldT = solnData(DELE_VAR,:,:,:)

     call Plasma_Solve(solnData(DELE_VAR,:,:,:), solnData(GNE_VAR,:,:,:), oldT, &
                       solnData(DFUN_VAR,:,:,:), solnData(DFEL_VAR,:,:,:),&
                       dt,del(DIR_X),del(DIR_Y),&
                       blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                       blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),T_res1)
     T_resBlock = T_resBlock + T_res1

     !obtain number density of heavy species
#ifdef DEBUG_PLASMA
     print *,"Going into Plasma solve heavy"
#endif
     do i=0,5
        oldT = solnData(DHV0_VAR+i,:,:,:)

        call Plasma_Solve(solnData(DHV0_VAR+i,:,:,:), solnData(GNH0_VAR+i,:,:,:), oldT, &
                          solnData(DFUN_VAR,:,:,:),   solnData(DFH0_VAR+i,:,:,:),&
                          dt,del(DIR_X),del(DIR_Y),&
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),T_res1)
        T_resBlockHV(i+1) = T_resBlockHV(i+1) + T_res1
     end do
     
     !obtain number density of ions
#ifdef DEBUG_PLASMA
     print *,"Going into Plasma solve ions"
#endif
     do i=0,3
        oldT = solnData(DHV6_VAR+i,:,:,:)
       
        call Plasma_Solve(solnData(DHV6_VAR+i,:,:,:),solnData(GNH6_VAR+i,:,:,:), oldT, &
                          solnData(DFUN_VAR,:,:,:),  solnData(DFEL_VAR,:,:,:),&
                          dt,del(DIR_X),del(DIR_Y),&
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),T_res1)
        T_resBlockHV(i+7) = T_resBlockHV(i+7) + T_res1
     end do

     k = 1
     do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
      do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)

           if(solnData(DFUN_VAR,i,j,k) .ge. 0.0) then

                solnData(DELE_VAR,i,j,k) = 0.99*1e18    ! Electrons
                solnData(DHV0_VAR,i,j,k) = 0.9*1e26     ! He
                solnData(DHV1_VAR,i,j,k) = 0.1*0.8*1e26 ! N2 
                solnData(DHV2_VAR,i,j,k) = 0.1*0.2*1e26 ! O2
                solnData(DHV6_VAR,i,j,k) = 0.9*1e18     ! He+
                solnData(DHV7_VAR,i,j,k) = 0.1*0.8*1e18 ! N2+
                solnData(DHV8_VAR,i,j,k) = 0.1*0.2*1e18 ! O2+
                solnData(DHV9_VAR,i,j,k) = 0.01*1e18    ! O-
                solnData(DNAT_VAR,i,j,k) = 1e26         !neutrals
                solnData(DNIT_VAR,i,j,k) = 1e18         !ions

           end if

       end do
     end do

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)


   end do

   if(blockCount .gt. 0) T_resBlock   = T_resBlock/blockCount
   if(blockCount .gt. 0) T_resBlockHV = T_resBlockHV/blockCount

   ! Collect residuals from other processes
   call MPI_Allreduce(T_resBlock, T_res, 1, FLASH_REAL,&
                      MPI_SUM, MPI_COMM_WORLD, ierr)
   
   call MPI_Allreduce(T_resBlockHV, T_resHV, 10, FLASH_REAL, &
                      MPI_SUM, MPI_COMM_WORLD, ierr)

   T_res   = T_res/pls_meshNumProcs
   T_resHV = T_resHV/pls_meshNumProcs

   if(pls_meshMe .eq. MASTER_PE) then
      
        print *,"Ne residual","",":",T_res

        do i=1,10
          print *,"Nh residual",i,":",T_resHV(i)
        end do

   endif

  gcMask = .FALSE.

  gcMask(DELE_VAR)  = .TRUE.
  gcMask(DFEL_VAR)  = .TRUE.
  gcMask(FVEA_VAR)  = .TRUE.
  gcMask(FVEI_VAR)  = .TRUE.
  gcMask(DNAT_VAR)  = .TRUE.
  gcMask(DNIT_VAR)  = .TRUE. 
  gcMask(TPHV_VAR)  = .TRUE.
  gcMask(TPEL_VAR)  = .TRUE.
  gcMask(GNE_VAR)   = .TRUE.
  gcMask(GNEBZ_VAR) = .TRUE.
  gcMask(GNERT_VAR) = .TRUE.

  do i=0,9
        gcMask(DHV0_VAR+i) = .TRUE.
        gcMask(DFH0_VAR+i) = .TRUE. 
        gcMask(GNH0_VAR+i) = .TRUE. 
  end do

  do i=0,16
        gcMask(RSP0_VAR+i) = .TRUE.
  end do

  call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)


end subroutine Plasma
