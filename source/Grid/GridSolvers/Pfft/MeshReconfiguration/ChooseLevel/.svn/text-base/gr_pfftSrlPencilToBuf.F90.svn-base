subroutine gr_pfftSrlPencilToBuf(item, buf)
#include "constants.h"
  use ut_conversionInterface, ONLY : ut_convertToMemoryOffset
  use gr_pfftFragmentNode, ONLY : fragment_node
  use gr_pfftData, ONLY : pfft_inLen
  use gr_pfftReconfigData, ONLY : pfft_pfftBuf
  implicit none
  type(fragment_node), pointer  :: item
  real, dimension(:), pointer :: buf
  integer :: index, i, j, k
  integer :: memoryOffset, oneDimensionalCoord
  integer, dimension(MDIM) :: elementCoord, arrayLBound, copyStart, copyEnd
  
  arrayLBound(:) = 1  !Required for ut_convertToMemoryOffset subroutine.
  copyStart(1:MDIM) = item % metadata % pfftStartPos
  copyEnd(1:MDIM) = item % metadata % pfftEndPos  
  index = item % metadata % bufStart

  !pfft_inLen is the equal amount of space assigned to each PFFT processor.
  do k = copyStart(KAXIS), copyEnd(KAXIS); elementCoord(KAXIS) = k
     do j = copyStart(JAXIS), copyEnd(JAXIS); elementCoord(JAXIS) = j
        do i = copyStart(IAXIS), copyEnd(IAXIS); elementCoord(IAXIS) = i
           call ut_convertToMemoryOffset(MDIM, elementCoord, arrayLBound, &
                pfft_inLen, memoryOffset)
           oneDimensionalCoord = memoryOffset + 1  !Fortran arrays 1-based.
           buf(index) = pfft_pfftBuf(oneDimensionalCoord)
           index = index + 1
        end do
     end do
  end do
end subroutine gr_pfftSrlPencilToBuf
