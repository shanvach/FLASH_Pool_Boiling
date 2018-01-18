subroutine Plasma( blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)

   use Plasma_data
   use Plasma_interface, only: Plasma_Solve, Plasma_hvDiffCoeff,&
                               Plasma_elDiffCoeff,Plasma_ColFreq

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
   integer :: ierr, i

   T_resBlock = 0.0

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
     
     !obtain diffusion coefficient of heavy particle species
     do i = 0,9
       call Plasma_hvDiffCoeff(solnData(DFH0_VAR+i,:,:,:), &
                               solnData(PRHV_VAR,:,:,:),solnData(TPHV_VAR,:,:,:), &
                               blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                               blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                               pls_RSCD(i+1), pls_MHSP(i+1), pls_MMIX(i+1)) 
     end do
     
     !obtain diffusion coefficient of electrons

     call Plasma_ColFreq(solnData(FVEI_VAR,:,:,:), solnData(FVEA_VAR,:,:,:), &
                         solnData(DELE_VAR,:,:,:), solnData(TPEL_VAR,:,:,:), &
                         blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                         blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS) )

     call Plasma_elDiffCoeff(solnData(DFEL_VAR,:,:,:),solnData(TPEL_VAR,:,:,:),&
                             solnData(FVEA_VAR,:,:,:),solnData(FVEI_VAR,:,:,:),&
                             blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                             blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS) )
     
     !obtain number density of electrons
     oldT = solnData(DELE_VAR,:,:,:)

     call Plasma_Solve(solnData(DELE_VAR,:,:,:), oldT, &
                       solnData(DFUN_VAR,:,:,:), solnData(DFEL_VAR,:,:,:),&
                       dt,del(DIR_X),del(DIR_Y),&
                       blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                       blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),T_res1)
     
     !obtain number density of heavy species
     do i=0,9
        oldT = solnData(DHV0_VAR+i,:,:,:)

        call Plasma_Solve(solnData(DHV0_VAR+i,:,:,:), oldT, &
                          solnData(DFUN_VAR,:,:,:), solnData(DFH0_VAR+i,:,:,:),&
                          dt,del(DIR_X),del(DIR_Y),&
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),T_res1)
     end do

     T_resBlock = T_resBlock + T_res1

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)


   end do

   if(blockCount .gt. 0) T_resBlock = T_resBlock/blockCount

   ! Collect residuals from other processes
   call MPI_Allreduce(T_resBlock, T_res, 1, FLASH_REAL,&
                      MPI_SUM, MPI_COMM_WORLD, ierr)

   T_res = T_res/pls_meshNumProcs

   if(pls_meshMe .eq. MASTER_PE) print *,"T_res:",T_res

  gcMask = .FALSE.

  gcMask(DELE_VAR)=.TRUE.
  gcMask(DFEL_VAR)=.TRUE.
  gcMask(FVEA_VAR)=.TRUE.
  gcMask(FVEI_VAR)=.TRUE.

  do i=0,9
        gcMask(DHV0_VAR+i) = .TRUE.
        gcMask(DFH0_VAR+i) = .TRUE.  
  end do

  call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)


end subroutine Plasma
