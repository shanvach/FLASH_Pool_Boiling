!!REORDER(4): solnData

subroutine gr_pfftSrlGridToBuf(item, buf)
#include "constants.h"
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
  use gr_pfftFragmentNode, ONLY : fragment_node
  use gr_pfftReconfigData, ONLY : pfft_srcGridVar
  implicit none
  type(fragment_node), pointer :: item
  real, dimension(:), pointer :: buf
  integer, dimension(MDIM) :: startCoords, endCoords
  real, dimension(:,:,:,:), pointer :: solnData
  integer :: gridVar, blockID, ii, jj, kk, index

  gridVar = pfft_srcGridVar
  blockID = item % metadata % srcFlashBlockID
  startCoords(1:MDIM) = item % metadata % srcFlashStartPos(1:MDIM)
  endCoords(1:MDIM) = item % metadata % srcFlashEndPos(1:MDIM)
  index = item % metadata % bufStart

  call Grid_getBlkPtr(blockID,solnData,CENTER)
  do kk = startCoords(KAXIS), endCoords(KAXIS)
     do jj = startCoords(JAXIS), endCoords(JAXIS)
        do ii = startCoords(IAXIS), endCoords(IAXIS)
           buf(index) = solnData(gridVar,ii,jj,kk)
           index = index + 1
        end do
     end do
  end do
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)
end subroutine gr_pfftSrlGridToBuf
