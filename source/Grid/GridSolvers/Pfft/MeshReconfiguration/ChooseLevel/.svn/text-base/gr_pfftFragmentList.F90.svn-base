!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/ChooseLevel/gr_pfftFragmentList
!!
!! NAME
!!  gr_pfftFragmentList
!!
!! SYNOPSIS
!!
!!  use gr_pfftFragmentList
!!   
!!***
module gr_pfftFragmentList
  use gr_pfftFragmentNode
  implicit none

  type fragment_list
     type(fragment_node), pointer :: H, T
  end type fragment_list

contains

  !An include file with the list operations.
#define CPP_NODE_NAME fragment_node
#define CPP_LIST_NAME fragment_list
#define CPP_NODE_DEFINITION \
  use gr_pfftFragmentNode, ONLY : fragment_node
#include "ut_listMethods.includeF90"

#undef CPP_NODE_NAME
#undef CPP_LIST_NAME
#undef CPP_NODE_DEFINITION

end module gr_pfftFragmentList
