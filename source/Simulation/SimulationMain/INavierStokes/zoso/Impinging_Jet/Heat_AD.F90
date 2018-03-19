subroutine Heat_AD( blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)

   use Heat_AD_interface, only: Heat_Solve
   use Grid_interface, only: Grid_getDeltas, Grid_getBlkIndexLimits,&
                             Grid_getBlkPtr, Grid_releaseBlkPtr,    &
                             Grid_fillGuardCells

   use IncompNS_data, only: ins_invRe, ins_meshMe, ins_meshNumProcs 

   implicit none
#include "constants.h"
#include "Heat_AD.h"
#include "Flash.h"   

   include "Flash_mpi.h"

   integer, INTENT(INOUT) :: blockCount
   integer, INTENT(IN) :: sweepOrder
   integer, INTENT(INOUT) :: blockList(MAXBLOCKS)
   real,    INTENT(IN) :: timeEndAdv, dt, dtOld

   integer ::  blockID,lb,ierr
   real ::  del(MDIM)
   integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
   logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)
   real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData
   real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: oldT

   real :: T_res, T_resblock, T_res1

   T_resblock = 0.0
   T_res      = 0.0
   T_res1     = 0.0

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

     oldT = solnData(TEMP_VAR,:,:,:)


     call Heat_Solve(solnData(TEMP_VAR,:,:,:), oldT,&
                     facexData(VELC_FACE_VAR,:,:,:),&
                     faceyData(VELC_FACE_VAR,:,:,:),&
                     facezData(VELC_FACE_VAR,:,:,:),&
                     dt,del(DIR_X),del(DIR_Y),del(DIR_Z),&
                     ins_invRe,&
                     blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                     blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                     blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),T_res1)

     T_resblock = T_resblock + T_res1

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

   end do

  if(blockCount .gt. 0) T_resBlock   = T_resBlock/blockCount

  ! Collect residuals from other processes
  call MPI_Allreduce(T_resBlock, T_res, 1, FLASH_REAL,&
                     MPI_SUM, MPI_COMM_WORLD, ierr)

  T_res = T_res/ins_meshNumProcs

  if(ins_meshMe .eq. MASTER_PE) print *,"Temperature residual: ",sqrt(T_res) 

  gcMask = .FALSE.
  gcMask(TEMP_VAR)=.TRUE.
  call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)


end subroutine Heat_AD
