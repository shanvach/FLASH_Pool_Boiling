!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/ChooseLevel/gr_pfftMapToInput
!!
!! NAME 
!!
!! gr_pfftMapToInput
!!
!! SYNOPSIS
!!
!! gr_pfftMapToInput(integer (IN) :: gridVar, &
!!                   real, target :: pfftInputArray(:))
!!
!! DESCRIPTION 
!!
!! Invokes lower level routines which copy data from the grid at variable 
!! "gridVar" and store it into the pencil input array.  This is non-trivial 
!! because the data is obtained from overlapping FLASH grid elements which 
!! may exist on different processors.
!!
!! ARGUMENTS
!!
!! gridVar - An index corresponding to the source grid variable.
!! pfftInputArray - The input pencil array.
!!
!! 
!!***

subroutine gr_pfftMapToInput(gridVar, pfftInputArray)
#include "constants.h"
#include "Flash.h"
#include "Pfft.h"
  use gr_pfftData, ONLY : pfft_comm
  use gr_pfftReconfigData, ONLY : pfft_pfftBuf, pfft_srcGridVar, &
       pfft_flashMsgList, pfft_pencilMsgList, pfft_mpiTag4, &
       pfft_logMode, pfft_logUnit
  use gr_pfftReconfigFn, ONLY : gr_pfftFnArgSendData, &
       gr_pfftFnArgRecvData, gr_pfftFnArgDeliverData
  use gr_pfftCommunicationPattern, ONLY : ut_mpiOverlapSendRecv, &
       common_state, send_state, recv_state
  use gr_pfftMessageList, ONLY : list_size
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  integer, intent(IN) :: gridVar
  real, dimension(:), target :: pfftInputArray
  include "Flash_mpi.h"

  real :: secondsTillAbort
  integer :: err
  logical :: useBarriers
  type(common_state), pointer :: commonState
  type(send_state), pointer :: sendState
  type(recv_state), pointer :: recvState

  pfft_pfftBuf => pfftInputArray
  pfft_srcGridVar = gridVar
  secondsTillAbort = 30.0
  useBarriers = .true.

  allocate(commonState, sendState, recvState, STAT=err)
  if (err /= 0) then
     call Driver_abortFlash("[gr_pfftMapToInput]: State allocation failed.")
  end if

  commonState % mpiTag = pfft_mpiTag4 
  commonState % mpiType = FLASH_REAL
  commonState % mpiComm = pfft_comm(IAXIS)
  commonState % direction = TO_PFFT

  sendState % msgList => pfft_flashMsgList
  nullify (sendState % activeMsg)
  sendState % numMsgSent = 0

  recvState % msgList => pfft_pencilMsgList
  nullify (recvState % activeMsg)
  recvState % numMsgRecv = 0
  recvState % allRecvPosted = .false.
  call list_size(pfft_pencilMsgList, recvState % totMsgToRecv)


  if (pfft_logMode .eqv. .true.) then
     write(pfft_logUnit,*) ""
     write(pfft_logUnit,*) "[gr_pfftMapToInput]: Communicating data"
     write(pfft_logUnit,*) "-------------------------------------------------------------"
  end if

  call ut_mpiOverlapSendRecv &
       (gr_pfftFnArgSendData, gr_pfftFnArgRecvData, &
       gr_pfftFnArgDeliverData, commonState, sendState, recvState, &
       secondsTillAbort, useBarriers)


  nullify (sendState % activeMsg, sendState % msgList, recvState % msgList)
  deallocate(commonState, sendState, recvState, STAT=err)
  if (err /= 0) then
     call Driver_abortFlash("[gr_pfftMapToInput]: State deallocation failed.")
  end if
  nullify(pfft_pfftBuf)

end subroutine gr_pfftMapToInput
