subroutine Heat_AD_init(blockCount,blockList)

   use Heat_AD_data
   use Grid_interface, only: Grid_getBlkPtr, Grid_releaseBlkPtr
   use IncompNS_data, only: ins_Ra, ins_Pr
   use RuntimeParameters_interface, only: RuntimeParameters_get

   implicit none

#include "constants.h"
#include "Flash.h"

   integer, INTENT(INOUT) :: blockCount
   integer, INTENT(INOUT) :: blockList(MAXBLOCKS)

   integer ::  blockID,lb
   real, pointer, dimension(:,:,:,:) :: solnData
   real :: rnd

   call RuntimeParameters_get('Twall_high', ht_Twall_high)
   call RuntimeParameters_get('Twall_low', ht_Twall_low)

   ht_invsqrtRaPr = 1. / (ins_Ra * ins_Pr)**0.5

   do lb = 1,blockCount
     blockID = blockList(lb)

     call random_number(rnd)

     call Grid_getBlkPtr(blockID,solnData,CENTER)
     solnData(TEMP_VAR,:,:,:) = rnd
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)

   end do

end subroutine Heat_AD_init
