subroutine gr_pfftFnArgSendMetadata(commonState, sendState, sendMore)
  use gr_pfftCommunicationState, ONLY : common_state, send_state
  use gr_pfftFragmentNode, ONLY : fragment_node, print_node
  use Driver_interface, ONLY : Driver_abortFlash, Driver_checkMPIErrorCode
  use gr_pfftData, ONLY : pfft_myPE
  use gr_pfftReconfigData, ONLY : pfft_logMode, pfft_logUnit
#include "constants.h"
#include "Flash.h"
#include "Pfft.h"

  implicit none
  include "Flash_mpi.h"
  type(common_state), pointer :: commonState
  type(send_state), pointer :: sendState
  logical, intent(OUT) :: sendMore

  type(fragment_node), pointer :: fragment
  integer :: destPE, ierr

  !sendState % activeMsg keeps track of our current message node in 
  !the send message data structure.
  if (sendState % numMsgSent == 0) then
     sendState % activeMsg => sendState % msgList % H
  else
     sendState % activeMsg => sendState % activeMsg % next
  end if  


  if (.not.(associated(sendState % activeMsg))) then
     sendMore = .false.
  else
     !A message object consists of multiple fragments to be sent to the 
     !same process.  We take the easy option of sending one fragment at 
     !a time.  This saves the effort of serialising the fragments into 
     !contiguous memory space.  Profiles should reveal that this easy 
     !option is not that bad since communication is non-blocking and 
     !the fragment messages are very small.

     fragment => sendState % activeMsg % fragmentList % H

     do while (associated(fragment))
        
        if (commonState % direction == TO_PFFT) then
           destPE = fragment % metadata % pfftProcID
        else if (commonState % direction == FROM_PFFT) then
           destPE = fragment % metadata % solnFlashProcID
        else
           call Driver_abortFlash("[gr_pfftFnArgSendMetadata]: "//&
                "Pfft directional movement error")
        end if

        !It is possible for process P to send to process P.  This is fine.
        call MPI_Isend(fragment % metadata, 1, commonState % mpiType, &
             destPE, commonState % mpiTag, commonState % mpiComm, &
             fragment % metadata_request, ierr)
        call Driver_checkMPIErrorCode(ierr)


        if (pfft_logMode .eqv. .true.) then
           write(pfft_logUnit,*) "[gr_pfftFnArgSendMetadata]: "//&
                "Sent following fragment to process ", destPE
           call print_node(pfft_logUnit, fragment)
        end if

        sendState % numMsgSent = sendState % numMsgSent + 1
        sendState % activeMsg % bufSize = &
             sendState % activeMsg % bufSize + &
             fragment % metadata % bufSize
        fragment => fragment % next
     end do

     sendMore = associated(sendState % activeMsg % next) 
  end if

end subroutine gr_pfftFnArgSendMetadata
