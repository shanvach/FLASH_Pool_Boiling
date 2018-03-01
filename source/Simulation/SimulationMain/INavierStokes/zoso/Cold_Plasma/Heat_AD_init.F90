subroutine Heat_AD_init(blockCount,blockList)

   use Heat_AD_data
   use Grid_interface, only: Grid_getBlkPtr, Grid_releaseBlkPtr
   use IncompNS_data, only: ins_invRe, ins_meshMe
   use RuntimeParameters_interface, ONLY : RuntimeParameters_get

   implicit none

#include "constants.h"
#include "Flash.h"

   include "Flash_mpi.h"

   integer, INTENT(INOUT) :: blockCount
   integer, INTENT(INOUT) :: blockList(MAXBLOCKS)

   integer ::  blockID,lb
   real, pointer, dimension(:,:,:,:) :: solnData

   call RuntimeParameters_get("Pr",ht_Pr)

   if (ins_meshMe .eq. MASTER_PE) then
     write(*,*) 'ht_Pr   =',ht_Pr
   end if

   do lb = 1,blockCount
     blockID = blockList(lb)

     call Grid_getBlkPtr(blockID,solnData,CENTER)
     solnData(TEMP_VAR,:,:,:) = 0.0
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)

   end do

end subroutine Heat_AD_init
