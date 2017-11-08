subroutine Heat_AD( blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder)

   use Heat_AD_interface, only: Heat_Solve
   use Grid_interface, only: Grid_getDeltas, Grid_getBlkIndexLimits,&
                             Grid_getBlkPtr, Grid_releaseBlkPtr,    &
                             Grid_fillGuardCells

   use IncompNS_data, only: ins_invRe,ins_alfa

   use ImBound_data, only: ib_temp_flg,ib_vel_flg

   use ImBound_interface, ONLY : ImBound

   implicit none
#include "constants.h"
#include "Heat_AD.h"
#include "Flash.h"   
#include "ImBound.h"

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

     oldT = solnData(TEMP_VAR,:,:,:)


     call Heat_Solve(solnData(TEMP_VAR,:,:,:), oldT,&
                     facexData(VELC_FACE_VAR,:,:,:),&
                     faceyData(VELC_FACE_VAR,:,:,:),&
                     dt,del(DIR_X),del(DIR_Y),del(DIR_Z),&
                     ins_invRe,&
                     blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                     blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS))

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)


   end do

  gcMask = .FALSE.
  gcMask(TEMP_VAR)=.TRUE.
  call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)

  ib_temp_flg = .true.
  ib_vel_flg  = .false.

  call ImBound( blockCount, blockList, ins_alfa*dt,FORCE_FLOW)

  gcMask = .FALSE.
  gcMask(TEMP_VAR)=.TRUE.
  call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)

end subroutine Heat_AD
