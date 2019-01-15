subroutine Heat_AD_init(blockCount,blockList)

   use Heat_AD_data
   use Grid_interface, only: Grid_getBlkPtr, Grid_releaseBlkPtr
   use IncompNS_data, only: ins_invRe
   use RuntimeParameters_interface, only: RuntimeParameters_get

   implicit none

#include "constants.h"
#include "Flash.h"

   integer, INTENT(INOUT) :: blockCount
   integer, INTENT(INOUT) :: blockList(MAXBLOCKS)

   integer ::  blockID,lb
   real, pointer, dimension(:,:,:,:) :: solnData

   ht_Pr = 0.7
   ht_Nu = 0.332*(ht_Pr**0.33)/(ins_invRe**0.5)
   ht_Bi = 1.0

   call RuntimeParameters_get('Twall_high', ht_Twall_high)
   call RuntimeParameters_get('Twall_low', ht_Twall_low)

   do lb = 1,blockCount
     blockID = blockList(lb)

     call Grid_getBlkPtr(blockID,solnData,CENTER)
     solnData(TEMP_VAR,:,:,:) = 0.0
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)

   end do

end subroutine Heat_AD_init
