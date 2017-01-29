!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/ChooseLevel/gr_pfftCommunicationState
!!
!! NAME
!!  gr_pfftCommunicationState
!!
!! SYNOPSIS
!!
!!  gr_pfftCommunicationState
!!   
!!***
module gr_pfftCommunicationState
  use gr_pfftMessageList, ONLY : message_list, message_node
  implicit none

  !Shared state: Read-only variables.
  type common_state     
     integer :: mpiTag, mpiType, mpiComm, direction
  end type common_state


  !Private state:
  type send_state
     type(message_list), pointer :: msgList
     type(message_node), pointer :: activeMsg  !Traverses message nodes.
     integer :: numMsgSent
  end type send_state

  type recv_state
     type(message_list), pointer :: msgList
     type(message_node), pointer :: activeMsg  !Traverses message nodes.
     integer :: numMsgRecv

     !For metadata exchange only:
     integer :: numGridPointsRecv, totGridPoints

     !For data exchange only:
     logical :: allRecvPosted
     integer :: totMsgToRecv
  end type recv_state

end module gr_pfftCommunicationState
