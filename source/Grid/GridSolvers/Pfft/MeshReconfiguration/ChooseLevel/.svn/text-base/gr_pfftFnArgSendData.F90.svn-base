subroutine gr_pfftFnArgSendData(commonState, sendState, sendMore)
  use gr_pfftCommunicationState, ONLY : common_state, send_state
  use Driver_interface, ONLY : Driver_abortFlash, Driver_checkMPIErrorCode
  use gr_pfftData, ONLY : pfft_myPE
  use gr_pfftReconfigData, ONLY : pfft_logMode, pfft_logUnit
  use gr_pfftReconfigFn, ONLY : gr_pfftDataToBuf
  use gr_pfftMessageNode, ONLY : print_node
#include "constants.h"
#include "Flash.h"
#include "Pfft.h"

  implicit none
  include "Flash_mpi.h"
  type(common_state), pointer :: commonState
  type(send_state), pointer :: sendState
  logical, intent(OUT) :: sendMore
  integer :: destPE, ierr


  if (sendState % numMsgSent == 0) then
     !This is executed the first time we call this subroutine.
     sendState % activeMsg => sendState % msgList % H
  end if  


  if (.not.(associated(sendState % activeMsg))) then
     sendMore = .false.  !Terminate if 0 messages in list.
  else

     !Each time we call this subroutine we only send 1 message, where the 
     !data for this 1 message can come from several different block fragments.
     if (pfft_logMode .eqv. .true.) then
        write(pfft_logUnit,*) "[gr_pfftFnArgSendData]: "//&
             "Sending following message to ", &
             sendState % activeMsg % PE_partner
        call print_node(pfft_logUnit, sendState % activeMsg)
     end if

     call gr_pfftDataToBuf(sendState % activeMsg, commonState % direction)

     call MPI_Isend(sendState % activeMsg % buf(1), &
          sendState % activeMsg % bufSize, &
          commonState % mpiType, &
          sendState % activeMsg % PE_partner, &
          commonState % mpiTag, &
          commonState % mpiComm, &
          sendState % activeMsg % data_request, &
          ierr)
     call Driver_checkMPIErrorCode(ierr)

     sendState % numMsgSent = sendState % numMsgSent + 1
     sendState % activeMsg => sendState % activeMsg % next
     sendMore = associated(sendState % activeMsg)
  end if

end subroutine gr_pfftFnArgSendData
