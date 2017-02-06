subroutine Heat_AD_init(blockCount,blockList)

   use Heat_AD_data
   use Grid_interface, only: Grid_getBlkPtr, Grid_releaseBlkPtr
   use IncompNS_data, only: ins_invRe,ins_meshMe
   use Multiphase_data, only: mph_thco2,mph_vis2,mph_cp2,mph_rho2,&
                              mph_thco1,mph_vis1,mph_cp1,mph_rho1

   use Grid_interface, ONLY : Grid_getDeltas,         &
                                 Grid_getBlkIndexLimits, &
                                 Grid_getCellCoords,     &
                                 Grid_getBlkPtr,         &
                                 Grid_releaseBlkPtr,     &
                                 Grid_getBlkBoundBox,    &
                                 Grid_getBlkCenterCoords

   use RuntimeParameters_interface, ONLY : RuntimeParameters_get

   implicit none

#include "constants.h"
#include "Flash.h"

   include "Flash_mpi.h"

   integer, INTENT(INOUT) :: blockCount
   integer, INTENT(INOUT) :: blockList(MAXBLOCKS)

   integer ::  blockID,lb,i,j,k
   real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData
   integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
   integer :: ierr,iter
   real :: maxdfun_local, maxdfun_global
   real :: beta, chi, soln, a_I, b_I, x1, x2, f1, f2, h


   call RuntimeParameters_get("Pr",ht_Pr)
   call RuntimeParameters_get("St",ht_St)


   if (ins_meshMe .eq. MASTER_PE) then

     write(*,*) 'ht_Pr   =',ht_Pr
     write(*,*) 'ht_St   =',ht_St

   end if

   ht_Twall_low  = 1.0
   ht_Twall_high = 0.0
   ht_Tsat       = 0.0

end subroutine Heat_AD_init
