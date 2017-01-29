subroutine gr_pfftFnArgDeliverMetadata(commonState, sendState, recvState)
#include "Pfft.h"
  use gr_pfftCommunicationState, ONLY : common_state, send_state, recv_state
  use gr_pfftMessageNode, ONLY : message_node
  use gr_pfftFragmentNode, ONLY : fragment_node
  use Driver_interface, ONLY : Driver_abortFlash, Driver_checkMPIErrorCode
  use gr_pfftReconfigData, ONLY : pfft_logMode, pfft_logUnit
  use gr_pfftData, ONLY : pfft_myPE
  implicit none  
  include "Flash_mpi.h"

  type(common_state), pointer :: commonState
  type(send_state), pointer :: sendState
  type(recv_state), pointer :: recvState

  type(message_node), pointer :: messageNode
  type(fragment_node), pointer :: fragmentNode
  integer :: destPE, ierr

  nullify(messageNode, fragmentNode)
  if ( (sendState % numMsgSent < 0) .or. (recvState % numMsgRecv < 0) ) then
     call Driver_abortFlash &
          ("[gr_pfftFnArgDeliverAllMeta]: Unexpected counter value")
  end if
  
  !We have already ensured delivery of all our received messages.  
  !We just need to ensure delivery of our sent messages.
  if (sendState % numMsgSent > 0) then

     messageNode => sendState % msgList % H
     do while(associated(messageNode))
        
        fragmentNode => messageNode % fragmentList % H
        do while(associated(fragmentNode))
           
           call MPI_Wait(fragmentNode % metadata_request, &
                fragmentNode % metadata_status, ierr)
           call Driver_checkMPIErrorCode(ierr)

           
           if (pfft_logMode .eqv. .true.) then
              if (commonState % direction == TO_PFFT) then
                 destPE = fragmentNode % metadata % pfftProcID
              else if (commonState % direction == FROM_PFFT) then
                 destPE = fragmentNode % metadata % solnFlashProcID
              else
                 call Driver_abortFlash("[gr_pfftFnArgSendMetadata]: "//&
                      "Pfft directional movement error")
              end if

              write(pfft_logUnit,*) "[gr_pfftFnArgDeliverMetadata]: "//&
                   "Delivered fragment to process ", destPE
           end if


           fragmentNode => fragmentNode % next
        end do

        messageNode => messageNode % next
     end do
  end if
  nullify(messageNode, fragmentNode)

end subroutine gr_pfftFnArgDeliverMetadata
