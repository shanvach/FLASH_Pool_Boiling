subroutine Heat_AD_init(blockCount,blockList)

   use Heat_AD_data

   use Grid_interface, only: Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkCenterCoords,&
                             Grid_getBlkIndexLimits

   use IncompNS_data, only: ins_Ra, ins_Pr, ins_isCoupled

   use RuntimeParameters_interface, only: RuntimeParameters_get

   implicit none

#include "constants.h"
#include "Flash.h"

   integer, INTENT(INOUT) :: blockCount
   integer, INTENT(INOUT) :: blockList(MAXBLOCKS)

   integer ::  blockID,lb
   integer, dimension(2,MDIM) :: blkLimits,blkLimitsGC 
   real, pointer, dimension(:,:,:,:) :: solnData
   real, dimension(MDIM) :: coord

   call RuntimeParameters_get('Twall_high', ht_Twall_high)
   call RuntimeParameters_get('Twall_low', ht_Twall_low)

   ht_invsqrtRaPr = 1. / (ins_Ra * ins_Pr)**0.5

   do lb = 1,blockCount
     blockID = blockList(lb)

     call Grid_getBlkCenterCoords(blockID,coord)     
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

#if NDIM == 3

     solnData(TEMP_VAR, blkLimitsGC(LOW,IAXIS)   :blkLimitsGC(HIGH,IAXIS),&
                        blkLimitsGC(LOW,JAXIS)   :blkLimitsGC(HIGH,JAXIS),&
                        blkLimitsGC(LOW,KAXIS)   :blkLimitsGC(HIGH,KAXIS))   = 1.0

     solnData(TEMP_VAR, blkLimitsGC(LOW,IAXIS)   :blkLimitsGC(HIGH,IAXIS),&
                        blkLimitsGC(LOW,JAXIS)   :blkLimitsGC(HIGH,JAXIS),&
                        blkLimitsGC(HIGH,KAXIS)/2:blkLimitsGC(HIGH,KAXIS)  ) = 0.0
#elif NDIM == 2
     IF (coord(JAXIS) <= 0.5) THEN
         solnData(TEMP_VAR,:,:,:) = 1.0
     ELSE
        solnData(TEMP_VAR,:,:,:) = 0.0
     END IF
#endif

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)

   end do

end subroutine Heat_AD_init
