subroutine Heat_AD( blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)

   use Heat_AD_interface, only: Heat_Solve2d, Heat_Solve3d
   use Grid_interface, only: Grid_getDeltas, Grid_getBlkIndexLimits,&
                             Grid_getBlkPtr, Grid_releaseBlkPtr,    &
                             Grid_fillGuardCells

   implicit none
#include "constants.h"
#include "Heat_AD.h"
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
#endif

     oldT = solnData(TEMP_VAR,:,:,:)

#if NDIM == 3
     call Heat_Solve3d(solnData(TEMP_VAR,:,:,:), oldT,&
                       facexData(VELC_FACE_VAR,:,:,:),&
                       faceyData(VELC_FACE_VAR,:,:,:),&
                       favezData(VELC_FACE_VAR,:,:,:),&
                       dt,del(DIR_X),del(DIR_Y),del(DIR_Z),&
                       blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                       blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                       blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS))
#elif NDIM == 2
     call Heat_Solve2d(solnData(TEMP_VAR,:,:,:), oldT,&
                       facexData(VELC_FACE_VAR,:,:,:),&
                       faceyData(VELC_FACE_VAR,:,:,:),&
                       dt,del(DIR_X),del(DIR_Y),del(DIR_Z),&
                       blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                       blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS))
#endif

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)

#if NDIM == 3
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

   end do

  gcMask = .FALSE.
  gcMask(TEMP_VAR)=.TRUE.
  call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)


end subroutine Heat_AD
