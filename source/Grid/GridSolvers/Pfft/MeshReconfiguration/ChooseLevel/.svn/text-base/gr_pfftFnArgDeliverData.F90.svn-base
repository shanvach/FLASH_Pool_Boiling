subroutine gr_pfftFnArgDeliverData(commonState, sendState, recvState)
#include "constants.h"
#include "Pfft.h"
  use gr_pfftCommunicationState, ONLY : common_state, send_state, recv_state
  use gr_pfftMessageNode, ONLY : message_node, print_node
  use Driver_interface, ONLY : Driver_abortFlash, Driver_checkMPIErrorCode
  use gr_pfftReconfigData, ONLY : pfft_logMode, pfft_logUnit
  use gr_pfftData, ONLY : pfft_myPE, pfft_comm
  implicit none  
  include "Flash_mpi.h"

  type(common_state), pointer :: commonState
  type(send_state), pointer :: sendState
  type(recv_state), pointer :: recvState

  type(message_node), pointer :: messageNode
  integer :: ierr

  nullify(messageNode)
  if ( (sendState % numMsgSent < 0) .or. (recvState % numMsgRecv < 0) ) then
     call Driver_abortFlash &
          ("[gr_pfftFnArgDeliverData]: ERROR - Unexpected counter value")
  end if
  
  !We have already ensured delivery of all our received messages.  
  !We just need to ensure delivery of our sent messages.
  if (sendState % numMsgSent > 0) then
     messageNode => sendState % msgList % H
     do while(associated(messageNode))   

        call MPI_Wait(messageNode % data_request, messageNode % data_status, ierr)
        call Driver_checkMPIErrorCode(ierr)

        if (pfft_logMode) then
           write(pfft_logUnit,*) "[gr_pfftFnArgDeliverData]: "//&
                "Delivered message to process ", messageNode % PE_partner
        end if


        !! This is important !!
        !! We now set messageNode % msgDelivered to .false. so that 
        !! it is in the correct state for the next time we 
        !! need to exchange data.
        messageNode % msgDelivered = .false.


        messageNode => messageNode % next
     end do

  end if
  nullify(messageNode)

end subroutine gr_pfftFnArgDeliverData
