subroutine gr_pfftFnArgRecvData(commonState, recvState, recvMore)
  use gr_pfftCommunicationState, ONLY : common_state, recv_state
  use gr_pfftMessageNode, ONLY : message_node, print_node
  use gr_pfftData, ONLY : pfft_myPE
  use Driver_interface, ONLY : Driver_checkMPIErrorCode, Driver_abortFlash
  use gr_pfftReconfigData, ONLY : pfft_logMode, pfft_logUnit
  use gr_pfftReconfigFn, ONLY : gr_pfftDataFromBuf
#include "constants.h"
#include "Flash.h"

  implicit none
  include "Flash_mpi.h"
  type(common_state), pointer :: commonState
  type(recv_state), pointer :: recvState
  logical, intent(OUT) :: recvMore

  type(message_node), pointer :: messageNode
  integer :: ierr
  logical :: msgNotifier


  !On the first call we post all receives.
  if (recvState % allRecvPosted .eqv. .false.) then
     messageNode => recvState % msgList % H
     do while (associated(messageNode))
        
        call MPI_Irecv(messageNode % buf(1), &
             messageNode % bufSize, &
             commonState % mpiType, &
             messageNode % PE_partner, &
             commonState % mpiTag, &
             commonState % mpiComm, &
             messageNode % data_request, &
             ierr)
        call Driver_checkMPIErrorCode(ierr)

        if (pfft_logMode) then
           write(pfft_logUnit,*) "[gr_pfftFnArgRecvData]: "//&
                "Preparing to receive from process ", messageNode % PE_partner
        end if

        messageNode => messageNode % next
     end do

     nullify(messageNode)  !We do not need messageNode anymore.
     recvState % allRecvPosted = .true.
     recvState % activeMsg => recvState % msgList % H
  end if


  !This loop cycles through the messages and if it finds a message 
  !that has been delivered it begins unpacking this data.  We only 
  !unpack a maximum of 1 message so that we distribute our time 
  !evenly between sending and receiving messages.
  msgNotifier = .false.
  do while (associated(recvState % activeMsg) .and. &
       (msgNotifier .eqv. .false.))
     
     if (recvState % activeMsg % msgDelivered .eqv. .false.) then        
        !This usage of MPI_Test is like doing an MPI_Testany on 
        !an array of requests, but instead on a list of requests.  If no 
        !messages have been delivered we leave this subroutine 
        !and end up in the send subroutine -- we could not do this with 
        !blocking MPI_Wait calls.  Also posting MPI_Tests keeps 
        !communication progressing (see paper).
        call MPI_Test(recvState % activeMsg % data_request, &
             msgNotifier, recvState % activeMsg % data_status, ierr)
        call Driver_checkMPIErrorCode(ierr)
     
        if (msgNotifier .eqv. .true.) then
           recvState % activeMsg % msgDelivered = .true.
           recvState % numMsgRecv = recvState % numMsgRecv + 1

           if (pfft_logMode) then
              write(pfft_logUnit,*) "[gr_pfftFnArgRecvData]: "//&
                   "Received following message from process ", &
                   recvState % activeMsg % PE_partner
              call print_node(pfft_logUnit, recvState % activeMsg)
           end if

           call gr_pfftDataFromBuf(recvState % activeMsg, commonState % direction)
        end if
     end if

     recvState % activeMsg => recvState % activeMsg % next
  end do
  


  recvMore = (recvState % numMsgRecv < recvState % totMsgToRecv)

  if (recvMore .eqv. .true.) then
     if (.not.(associated(recvState % activeMsg))) then
        !This will happen if we have tested each node but not all messages 
        !have been delivered into all nodes.  So point to the head of the 
        !list and start the MPI_Test cycle again.
        recvState % activeMsg => recvState % msgList % H
     end if
  end if

end subroutine gr_pfftFnArgRecvData
