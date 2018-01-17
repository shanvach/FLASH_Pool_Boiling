subroutine Plasma( blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)

   use Plasma_data
   use Plasma_interface, only: Plasma_Solve
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
   integer :: ierr

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

     oldT = solnData(DELE_VAR,:,:,:)

     call Plasma_Solve(solnData(DELE_VAR,:,:,:), oldT, solnData(DFUN_VAR,:,:,:),&
                     dt,del(DIR_X),del(DIR_Y),&
                     blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                     blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),T_res1)

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
  call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)


end subroutine Plasma
