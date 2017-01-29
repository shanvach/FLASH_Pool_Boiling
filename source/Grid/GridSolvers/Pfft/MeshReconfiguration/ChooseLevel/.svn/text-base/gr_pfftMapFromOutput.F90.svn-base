!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/ChooseLevel/gr_pfftMapFromOutput
!!
!! NAME 
!!
!! gr_pfftMapFromOutput
!!
!! SYNOPSIS
!!
!! gr_pfftMapFromOutput(integer (IN) :: gridVar, &
!!                      real, target :: pfftOutputArray(:))
!!
!! DESCRIPTION 
!!
!! Invokes lower level routines which copy data from the pencil output array 
!! and store it in the grid at variable "gridVar".  This is non-trivial 
!! because the data is must be placed in overlapping FLASH grid elements which 
!! may exist on different processors.
!!
!! ARGUMENTS
!!
!! gridVar - An index corresponding to the source grid variable.
!! pfftOutputArray - The output pencil array.
!!
!!***

subroutine gr_pfftMapFromOutput(gridVar, pfftOutputArray)
#include "constants.h"
#include "Flash.h"
#include "Pfft.h"
  use gr_pfftData, ONLY : pfft_comm
  use gr_pfftReconfigData, ONLY : pfft_pfftBuf, pfft_solnGridVar, &
       pfft_flashMsgList, pfft_pencilMsgList, pfft_mpiTag5, pfft_logMode, &
       pfft_logUnit
  use gr_pfftCommunicationPattern, ONLY : ut_mpiOverlapSendRecv, &
       common_state, send_state, recv_state
  use gr_pfftReconfigFn, ONLY : gr_pfftFnArgSendData, &
       gr_pfftFnArgRecvData, gr_pfftFnArgDeliverData
  use gr_pfftMessageList, ONLY : list_size
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  integer, intent(IN) :: gridVar
  real, dimension(:), target :: pfftOutputArray
  include "Flash_mpi.h"

  real :: secondsTillAbort
  logical :: useBarriers
  type(common_state), pointer :: commonState
  type(send_state), pointer :: sendState
  type(recv_state), pointer :: recvState
  integer :: err

  pfft_pfftBuf => pfftOutputArray
  pfft_solnGridVar = gridVar
  secondsTillAbort = 30.0
  useBarriers = .true.

  allocate(commonState, sendState, recvState, STAT=err)
  if (err /= 0) then
     call Driver_abortFlash("[gr_pfftMapFromOutput]: State allocation failed.")
  end if
 
  commonState % mpiTag = pfft_mpiTag5
  commonState % mpiType = FLASH_REAL
  commonState % mpiComm = pfft_comm(IAXIS)
  commonState % direction = FROM_PFFT

  sendState % msgList => pfft_pencilMsgList
  nullify (sendState % activeMsg)
  sendState % numMsgSent = 0

  recvState % msgList => pfft_flashMsgList
  nullify (recvState % activeMsg)
  recvState % numMsgRecv = 0
  recvState % allRecvPosted = .false.
  call list_size(pfft_flashMsgList, recvState % totMsgToRecv)


  if (pfft_logMode .eqv. .true.) then
     write(pfft_logUnit,*) ""
     write(pfft_logUnit,*) "[gr_pfftMapFromOutput]: Communicating data"
     write(pfft_logUnit,*) "-------------------------------------------------------------"
  end if

  call ut_mpiOverlapSendRecv &
       (gr_pfftFnArgSendData, gr_pfftFnArgRecvData, &
       gr_pfftFnArgDeliverData, commonState, sendState, recvState, &
       secondsTillAbort, useBarriers)


  nullify (sendState % activeMsg, sendState % msgList, recvState % msgList)
  deallocate(commonState, sendState, recvState, STAT=err)
  if (err /= 0) then
     call Driver_abortFlash("[gr_pfftMapFromInput]: State deallocation failed.")
  end if
  nullify(pfft_pfftBuf)
end subroutine gr_pfftMapFromOutput
