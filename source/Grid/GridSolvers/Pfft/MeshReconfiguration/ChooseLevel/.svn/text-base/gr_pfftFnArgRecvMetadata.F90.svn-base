subroutine gr_pfftFnArgRecvMetadata(commonState, recvState, recvMore)
  use gr_pfftCommunicationState, ONLY : common_state, recv_state
  use gr_pfftMessageNode, ONLY : message_node
  use gr_pfftFragmentList, ONLY : create_node, push_back, print_node, &
       fragment_node
  use gr_pfftData, ONLY : pfft_myPE
  use Driver_interface, ONLY : Driver_checkMPIErrorCode, Driver_abortFlash
  use gr_pfftReconfigData, ONLY : pfft_pencilSize, pfft_logMode, pfft_logUnit
  use gr_pfftReconfigFn, ONLY : gr_pfftGetMessageNode

#include "constants.h"
#include "Flash.h"

  implicit none
  include "Flash_mpi.h"
  type(common_state), pointer :: commonState
  type(recv_state), pointer :: recvState
  logical, intent(OUT) :: recvMore
  integer, dimension(MPI_STATUS_SIZE) :: statusScalar

  type(message_node), pointer :: messageNode
  type(fragment_node), pointer :: fragmentNode
  integer :: pfftFragmentSize, ierr, srcMsgSize, srcProc
  logical :: msgNotifier

  nullify(messageNode, fragmentNode)
  recvMore = (recvState % totGridPoints > recvState % numGridPointsRecv)

  ProbeForMessages: if (recvMore .eqv. .true.) then

     !Message discovery
     !--------------------------------------------------------------------------
     !Probe for fragments.
     call MPI_Iprobe(MPI_ANY_SOURCE, commonState % mpiTag, &
          commonState % mpiComm, msgNotifier, statusScalar, ierr)
     call Driver_checkMPIErrorCode(ierr)
     !--------------------------------------------------------------------------

     !Message receive
     !--------------------------------------------------------------------------
     ReceiveMessage: if (msgNotifier .eqv. .true.) then
        srcProc = statusScalar(MPI_SOURCE)
        call MPI_Get_count(statusScalar, commonState % mpiType, srcMsgSize, ierr)
        call Driver_checkMPIErrorCode(ierr)
        if (srcMsgSize /= 1) then
           print *, "Processor:", pfft_myPE, "message tag:", &
                commonState % mpiTag, "message size:", srcMsgSize
           call Driver_abortFlash &
                ("[gr_pfftFnArgPossibleRecv]: Unexpected message size.")
        end if


        !Get a pointer to the message object which holds the individual metadata 
        !fragments from srcProc.  The returned message node is always valid.
        call gr_pfftGetMessageNode(srcProc, recvState % msgList, messageNode)

        call create_node(fragmentNode)
        call MPI_Recv(fragmentNode % metadata, 1, commonState % mpiType, &
             srcProc, commonState % mpiTag, commonState % mpiComm, &
             fragmentNode % metadata_status, ierr) !Must be blocking!
        call Driver_checkMPIErrorCode(ierr)

        call push_back(messageNode % fragmentList, fragmentNode)

        if (pfft_logMode .eqv. .true.) then
           write(pfft_logUnit,*) "[gr_pfftFnArgRecvMetadata]: "//&
                "Received following fragment from process ", srcProc
           call print_node(pfft_logUnit, fragmentNode)
           write(pfft_logUnit,"(a,i4,a,i10,a,i10)") &
                " [gr_pfftFnArgRecvMetadata]: Added fragment to message ", &
                messageNode % ID, " with data start ", fragmentNode % metadata % bufStart, &
                " and data size ", fragmentNode % metadata % bufSize
        end if

        messageNode % bufSize = messageNode % bufSize + &
             fragmentNode % metadata % bufSize

        !We accumulate the data message sizes into the "recvState % gridPoints" 
        !variable.  When this value reaches the size of the pencil grid 
        !assigned to myPE we know that we have received all metadata messages
        recvState % numGridPointsRecv = recvState % numGridPointsRecv + &
             fragmentNode % metadata % bufSize

        !Advance the receive message counter.
        recvState % numMsgRecv = recvState % numMsgRecv + 1

        !Note: Do not deallocate node memory because the metadata sits in this 
        !space!  We can retrieve this metadata later by following the head or 
        !tail list pointers.
        nullify(messageNode, fragmentNode)

     end if ReceiveMessage
     !--------------------------------------------------------------------------
  end if ProbeForMessages

end subroutine 
