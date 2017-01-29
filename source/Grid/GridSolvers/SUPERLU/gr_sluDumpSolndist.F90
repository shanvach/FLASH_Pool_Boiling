
#include "Flash.h"
#include "constants.h"

subroutine gr_sluDumpSolndist(isoln,nloc_dofs,solA)

  
  use Grid_interface,   ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr,   &
                               Grid_getBlkIndexLimits

  use superlu_common, only : neighProc_blkcnt,neighProc_blkList

  use Grid_data, only : gr_meshMe

  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  integer, intent(in) :: isoln,nloc_dofs
  real, intent(in)   :: solA(nloc_dofs)

  ! Local Variables
  integer :: blockList_set(MAXBLOCKS),blockCount_set
  integer,dimension(2,MDIM)::blkLimits,blkLimitsGC
  real, pointer, dimension(:,:,:,:) :: solnData

  integer :: i,ii,j,k,lb,blockID

  blockCount_set = neighProc_blkcnt(CONSTANT_TWO,CONSTANT_ONE)      ! subset blocks number in mype
  blockList_set(1:blockCount_set)=neighProc_blkList(1:blockCount_set,CONSTANT_ONE) ! subset blocks in mype

  ii = 0
  do lb=1,blockCount_set
     
     blockID = blockList_set(lb) ! This is in Processors

     ! Point to Blocks centered variables:
     call Grid_getBlkPtr(blockID,solnData,CENTER)

     ! Get Block limits
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              ii = ii + 1
              solnData(isoln,i,j,k) = solA(ii)
           enddo
        enddo
     enddo

     ! Release pointer
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  enddo

  if (ii .ne. nloc_dofs) &
  call Driver_abortFlash('gr_sluDumpSolndist : final ii .ne. nloc_dofs')

  return

end subroutine gr_sluDumpSolndist
