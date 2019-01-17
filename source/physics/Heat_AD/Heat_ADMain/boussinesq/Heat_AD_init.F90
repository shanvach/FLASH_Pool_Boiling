subroutine Heat_AD_init(blockCount,blockList)

   use Heat_AD_data
   use Grid_interface, only: Grid_getBlkPtr, Grid_releaseBlkPtr, &
                             Grid_getBlkCenterCoords
   use IncompNS_data, only: ins_Ra, ins_Pr
   use RuntimeParameters_interface, only: RuntimeParameters_get

   implicit none

#include "constants.h"
#include "Flash.h"

   integer, INTENT(INOUT) :: blockCount
   integer, INTENT(INOUT) :: blockList(MAXBLOCKS)

   integer ::  blockID,lb
   real, pointer, dimension(:,:,:,:) :: solnData
   real, dimension(MDIM) :: coord

   call RuntimeParameters_get('Twall_high', ht_Twall_high)
   call RuntimeParameters_get('Twall_low', ht_Twall_low)

   ht_invsqrtRaPr = 1. / (ins_Ra * ins_Pr)**0.5

   do lb = 1,blockCount
     blockID = blockList(lb)

     call Grid_getBlkCenterCoords(blockID,coord)     
     call Grid_getBlkPtr(blockID,solnData,CENTER)

#if NDIM == 3
     IF (coord(KAXIS) >= 0.5) THEN
         solnData(TEMP_VAR,:,:,:) = 0.0
     ELSE
        solnData(TEMP_VAR,:,:,:) = 1.0
     END IF
#elif NDIM == 2
     IF (coord(JAXIS) >= 0.5) THEN
         solnData(TEMP_VAR,:,:,:) = 0.0
     ELSE
        solnData(TEMP_VAR,:,:,:) = 1.0
     END IF
#endif

     call Grid_releaseBlkPtr(blockID,solnData,CENTER)

   end do

end subroutine Heat_AD_init
