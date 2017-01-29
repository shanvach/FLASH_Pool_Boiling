!Functions accept objects as arguments, so this module 
!only compiles when objects are defined.  Hence we cannot 
!place these functions in localAPI.

module gr_pfftReconfigFn

  interface
     subroutine gr_pfftFnArgSendMetadata(commonState, sendState, sendMore)
       use gr_pfftCommunicationState, ONLY : common_state, send_state
       implicit none
       type(common_state), pointer :: commonState
       type(send_state), pointer :: sendState
       logical, intent(OUT) :: sendMore
     end subroutine gr_pfftFnArgSendMetadata
  end interface

  interface
     subroutine gr_pfftFnArgRecvMetadata(commonState, recvState, recvMore)
       use gr_pfftCommunicationState, ONLY : common_state, recv_state
       implicit none
       type(common_state), pointer :: commonState
       type(recv_state), pointer :: recvState
       logical, intent(OUT) :: recvMore
     end subroutine gr_pfftFnArgRecvMetadata
  end interface

  interface
     subroutine gr_pfftFnArgDeliverMetadata(commonState, sendState, recvState)
       use gr_pfftCommunicationState, ONLY : common_state, send_state, recv_state
       implicit none  
       type(common_state), pointer :: commonState
       type(send_state), pointer :: sendState
       type(recv_state), pointer :: recvState
     end subroutine gr_pfftFnArgDeliverMetadata
  end interface

  interface
     subroutine gr_pfftFnArgSendData(commonState, sendState, sendMore)
       use gr_pfftCommunicationState, ONLY : common_state, send_state
       implicit none
       type(common_state), pointer :: commonState
       type(send_state), pointer :: sendState
       logical, intent(OUT) :: sendMore
     end subroutine gr_pfftFnArgSendData
  end interface

  interface
     subroutine gr_pfftFnArgRecvData(commonState, recvState, recvMore)
       use gr_pfftCommunicationState, ONLY : common_state, recv_state
       implicit none
       type(common_state), pointer :: commonState
       type(recv_state), pointer :: recvState
       logical, intent(OUT) :: recvMore
     end subroutine gr_pfftFnArgRecvData
  end interface

  interface
     subroutine gr_pfftFnArgDeliverData(commonState, sendState, recvState)
       use gr_pfftCommunicationState, ONLY : common_state, send_state, recv_state
       implicit none  
       type(common_state), pointer :: commonState
       type(send_state), pointer :: sendState
       type(recv_state), pointer :: recvState
     end subroutine gr_pfftFnArgDeliverData
  end interface

  interface
     subroutine gr_pfftDataToBuf(messageNode, direction)
       use gr_pfftMessageNode, ONLY : message_node
       implicit none
       type(message_node), pointer  :: messageNode
       integer, intent(IN) :: direction
     end subroutine gr_pfftDataToBuf
  end interface
  
  interface
     subroutine gr_pfftDataFromBuf(messageNode, direction)
       use gr_pfftMessageNode, ONLY : message_node
       implicit none
       type(message_node), pointer  :: messageNode
       integer, intent(IN) :: direction
     end subroutine gr_pfftDataFromBuf
  end interface
  
  interface
     subroutine gr_pfftSrlGridToBuf(item, buf)
       use gr_pfftFragmentNode, ONLY : fragment_node
       implicit none
       type(fragment_node), pointer :: item
       real, dimension(:), pointer :: buf
     end subroutine gr_pfftSrlGridToBuf
  end interface
  
  interface
     subroutine gr_pfftSrlPencilToBuf(item, buf)
       use gr_pfftFragmentNode, ONLY : fragment_node
       implicit none
       type(fragment_node), pointer  :: item
       real, dimension(:), pointer :: buf
     end subroutine gr_pfftSrlPencilToBuf
  end interface
  
  interface
     subroutine gr_pfftSrlBufToPencil(item, buf)
       use gr_pfftFragmentNode, ONLY : fragment_node
       implicit none
       type(fragment_node), pointer  :: item
       real, dimension(:), pointer :: buf
     end subroutine gr_pfftSrlBufToPencil
  end interface
  
  interface
     subroutine gr_pfftSrlBufToGrid(item, buf)
       use gr_pfftFragmentNode, ONLY : fragment_node
       implicit none
       type(fragment_node), pointer  :: item
       real, dimension(:), pointer :: buf
     end subroutine gr_pfftSrlBufToGrid
  end interface

  interface
     subroutine gr_pfftGetMessageNode(PE, messageList, messageNode)
       use gr_pfftMessageNode, ONLY : message_node
       use gr_pfftMessageList, ONLY : message_list
       implicit none
       integer, intent(IN) :: PE
       type(message_list), pointer :: messageList
       type(message_node), pointer :: messageNode
     end subroutine gr_pfftGetMessageNode
  end interface
end module gr_pfftReconfigFn
