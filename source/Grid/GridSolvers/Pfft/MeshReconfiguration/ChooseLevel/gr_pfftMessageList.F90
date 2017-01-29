!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/ChooseLevel/gr_pfftMessageList
!!
!! NAME
!!  gr_pfftMessageList
!!
!! SYNOPSIS
!!
!!  use gr_pfftMessageList
!!   
!!***
module gr_pfftMessageList
  use gr_pfftMessageNode
  implicit none

  type message_list
     type(message_node), pointer :: H, T
  end type message_list

contains

  !An include file with the list operations.
#define CPP_NODE_NAME message_node
#define CPP_LIST_NAME message_list
#define CPP_NODE_DEFINITION \
  use gr_pfftMessageNode, ONLY : message_node
#include "ut_listMethods.includeF90"

#define CPP_DATA_ARG \
  integer, intent(IN) :: dataArg
#include "ut_listFindNodeMethod.includeF90"

#undef CPP_NODE_NAME
#undef CPP_LIST_NAME
#undef CPP_NODE_DEFINITION
#undef CPP_DATA_ARG
end module gr_pfftMessageList
