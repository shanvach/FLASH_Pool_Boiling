!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/ChooseLevel/gr_pfftCommunicationPattern
!!
!! NAME
!!  gr_pfftCommunicationPattern
!!
!! SYNOPSIS
!!
!!  gr_pfftCommunicationPattern
!!   
!!***
module gr_pfftCommunicationPattern
  use gr_pfftCommunicationState, ONLY : send_state, recv_state, common_state
  implicit none
  include "Flash_mpi.h"

contains

#define CPP_COMMONSTATE_DEFINITION \
  use gr_pfftCommunicationState, ONLY : common_state
#define CPP_SENDSTATE_DEFINITION \
  use gr_pfftCommunicationState, ONLY : send_state
#define CPP_RECVSTATE_DEFINITION \
  use gr_pfftCommunicationState, ONLY : recv_state

  !An include file containing the communication pattern.
#include "ut_mpiOverlapSendRecv.includeF90"
#undef CPP_COMMONSTATE_DEFINITION
#undef CPP_SENDSTATE_DEFINITION
#undef CPP_RECVSTATE_DEFINITION

end module gr_pfftCommunicationPattern
