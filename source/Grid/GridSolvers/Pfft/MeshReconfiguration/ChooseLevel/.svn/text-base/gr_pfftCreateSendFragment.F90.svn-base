!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/ChooseLevel/gr_pfftCreateSendFragment
!!
!! NAME 
!!
!! gr_pfftCreateSendFragment
!!
!! SYNOPSIS
!!
!!
!! DESCRIPTION 
!!
!!
!! ARGUMENTS
!!
!!***

subroutine gr_pfftCreateSendFragment(flashProcID, flashBlockID, &
     lBlockFragStart, lBlockFragEnd, lActualBlockFragStart, lActualBlockFragEnd, &
     blkType, blkRefLev, solveLevel, pfftProcID, &
     lPencilFragStart, lPencilFragEnd) 

#include "constants.h"
#include "Flash.h"
  use gr_pfftMessageList, ONLY : message_node
  use gr_pfftFragmentList, ONLY : fragment_node, &
       create_fragment_node => create_node, push_back_fragment_node => push_back
  use gr_pfftReconfigData, ONLY : pfft_logMode, pfft_logUnit, pfft_flashMsgList
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_pfftReconfigFn, ONLY : gr_pfftGetMessageNode

  implicit none
  integer, intent(IN) :: flashProcID, flashBlockID
  integer, dimension(1:MDIM), intent(IN) :: lBlockFragStart, lBlockFragEnd, &
       lActualBlockFragStart, lActualBlockFragEnd
  integer, intent(IN) :: blkType, blkRefLev, solveLevel, pfftProcID
  integer, dimension(1:MDIM), intent(IN) :: lPencilFragStart, lPencilFragEnd 
  integer :: flashFragmentSize, pfftFragmentSize
  type(fragment_node), pointer :: fragmentNode
  type(message_node), pointer :: messageNode
  character(len=*), parameter  :: logStr = &
       "(A, I3, A, I3, A, I3, A, I3, A, I3, A, I3, A, I3, A, I3, A, I3, A, I3, A, I3, A, I3)"

  nullify(fragmentNode, messageNode)

  !Perform a consistency check:
  flashFragmentSize = &
       product( (lBlockFragEnd(1:MDIM) - lBlockFragStart(1:MDIM) + 1) )

  pfftFragmentSize = &
       product( (lPencilFragEnd(1:MDIM) - lPencilFragStart(1:MDIM) + 1) )

  if (flashFragmentSize /= pfftFragmentSize) then
     call Driver_abortFlash("[gr_pfftCreateSrlSendNode]: "//&
          "Mismatch in fragment calculation!")
  end if


  call create_fragment_node(fragmentNode)

  fragmentNode % metadata % srcFlashProcID = flashProcID
  fragmentNode % metadata % srcFlashBlockID = flashBlockID
  fragmentNode % metadata % srcFlashBlockType = blkType
  fragmentNode % metadata % srcFlashBlockRefLev = blkRefLev
  fragmentNode % metadata % solnFlashProcID = NONEXISTENT
  fragmentNode % metadata % solnFlashBlockID = NONEXISTENT
  fragmentNode % metadata % pfftProcID = pfftProcID

  fragmentNode % metadata % srcFlashStartPos = lBlockFragStart
  fragmentNode % metadata % srcFlashEndPos = lBlockFragEnd
  fragmentNode % metadata % srcFlashActualStartPos = lActualBlockFragStart
  fragmentNode % metadata % srcFlashActualEndPos = lActualBlockFragEnd
  fragmentNode % metadata % pfftStartPos = lPencilFragStart
  fragmentNode % metadata % pfftEndPos = lPencilFragEnd


  !We now have a node describing a block fragment. Figure
  !out whether to put it in a new message or in a pre-existing message.
  !gr_pfftGetMessageNode always returns a valid messageNode.
  call gr_pfftGetMessageNode(pfftProcID, pfft_flashMsgList, messageNode)
  call push_back_fragment_node(messageNode % fragmentList, fragmentNode)


  if (blkRefLev == solveLevel) then
     fragmentNode % metadata % bufSize = product(&
             fragmentNode % metadata % srcFlashEndPos(1:MDIM) - &
             fragmentNode % metadata % srcFlashStartPos(1:MDIM) + 1)
  else
     call Driver_abortFlash("[gr_pfftCreateSendFragment]: Not coded")
  end if


  !Multiple fragments' data will exist in the same message buffer.
  !We must ensure the bufStart/bufSize values do not overlap.
  if (associated(fragmentNode % prev)) then
     fragmentNode % metadata % bufStart = &
          fragmentNode % prev % metadata % bufStart + &
          fragmentNode % prev % metadata % bufSize
  end if


  if (pfft_logMode .eqv. .true.) then
     write(pfft_logUnit,*) "[gr_pfftCreateSendFragment]: "//&
          "Added fragment to message ", messageNode % ID, &
          "destined for process ", pfftProcID
  end if

  nullify(fragmentNode, messageNode)

end subroutine gr_pfftCreateSendFragment
