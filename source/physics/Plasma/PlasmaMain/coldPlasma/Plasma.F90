subroutine Plasma( blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)

   use Plasma_data
   use Plasma_interface, only: Plasma_Solve, Plasma_hvDiffCoeff,     &
                               Plasma_elDiffCoeff,Plasma_ColFreq,    &
                               Plasma_sumNeutrals,Plasma_spReactions,&
                               Plasma_spGeneration

   use Grid_interface, only: Grid_getDeltas, Grid_getBlkIndexLimits,&
                             Grid_getBlkPtr, Grid_releaseBlkPtr,    &
                             Grid_fillGuardCells
   implicit none
#include "constants.h"
#include "Plasma.h"
#include "Flash.h"   

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
   integer :: ierr, i

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
     
     !obtain total number of neutral particles
     ! important: We need to make sure that DHVT is set to zero at this point
     ! always, otherwise we will be mistakingly adding all neutral values from
     ! all time steps immediately below
     solnData(DHVT_VAR,:,:,:) = 0.0
     do i = 0,9
       call Plasma_sumNeutrals( solnData(DHV0_VAR+i,:,:,:),solnData(DHVT_VAR,:,:,:), &
                                blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                                blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS))
     end do 
     !obtain diffusion coefficient of heavy particle species
     do i = 0,9
       call Plasma_hvDiffCoeff(solnData(DFH0_VAR+i,:,:,:), &
                               solnData(PRHV_VAR,:,:,:),solnData(TPHV_VAR,:,:,:), &
                               blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                               blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                               pls_RSCD(i+1), pls_MHSP(i+1)) 
     end do
     
     !obtain collision frequencies
     call Plasma_ColFreq(solnData(FVEI_VAR,:,:,:), solnData(FVEA_VAR,:,:,:), &
                         solnData(DHVT_VAR,:,:,:), solnData(DELE_VAR,:,:,:), &
                         solnData(TPEL_VAR,:,:,:), &
                         blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                         blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS) )
     !obtain diffusion coefficient of electros
     call Plasma_elDiffCoeff(solnData(DFEL_VAR,:,:,:),solnData(TPEL_VAR,:,:,:),&
                             solnData(FVEA_VAR,:,:,:),solnData(FVEI_VAR,:,:,:),&
                             blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                             blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS) )
      
     !obtain chemical reaction rates for all species
     call Plasma_spReactions(solnData(RSP0_VAR,:,:,:), solnData(RSP1_VAR,:,:,:),&
                             solnData(RSP2_VAR,:,:,:), solnData(RSP3_VAR,:,:,:),&
                             solnData(RSP4_VAR,:,:,:), solnData(RSP5_VAR,:,:,:),&
                             solnData(RSP6_VAR,:,:,:), solnData(RSP7_VAR,:,:,:),&
                             solnData(RSP8_VAR,:,:,:), solnData(RSP9_VAR,:,:,:),&
                             solnData(RSP10_VAR,:,:,:), solnData(RSP11_VAR,:,:,:),&
                             solnData(RSP12_VAR,:,:,:), solnData(RSP13_VAR,:,:,:),&
                             solnData(TPHV_VAR,:,:,:), solnData(TPEL_VAR,:,:,:),&
                             blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                             blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS) )

     !obtain species generation rates using chemical reations
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
     oldT = solnData(DELE_VAR,:,:,:)

     call Plasma_Solve(solnData(DELE_VAR,:,:,:), solnData(GNE_VAR,:,:,:), oldT, &
                       solnData(DFUN_VAR,:,:,:), solnData(DFEL_VAR,:,:,:),&
                       dt,del(DIR_X),del(DIR_Y),&
                       blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                       blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),T_res1)
     
     T_resBlock = T_resBlock + T_res1

     !obtain number density of heavy species
     do i=0,9
        oldT = solnData(DHV0_VAR+i,:,:,:)

        call Plasma_Solve(solnData(DHV0_VAR+i,:,:,:), solnData(GNH0_VAR+i,:,:,:), oldT, &
                          solnData(DFUN_VAR,:,:,:),   solnData(DFH0_VAR+i,:,:,:),&
                          dt,del(DIR_X),del(DIR_Y),&
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),T_res1)
      
        T_resBlockHV(i+1) = T_resBlockHV(i+1) + T_res1

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

  gcMask(DELE_VAR)=.TRUE.
  gcMask(DFEL_VAR)=.TRUE.
  gcMask(FVEA_VAR)=.TRUE.
  gcMask(FVEI_VAR)=.TRUE.
  gcMask(DHVT_VAR)=.TRUE.

  do i=0,9
        gcMask(DHV0_VAR+i) = .TRUE.
        gcMask(DFH0_VAR+i) = .TRUE.  
  end do

  call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)


end subroutine Plasma
