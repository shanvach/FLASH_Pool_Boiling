!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/ChooseLevel/gr_pfftGenMap
!!
!! NAME
!!
!!  gr_pfftGenMap
!!
!! SYNOPSIS
!!  
!!  gr_pfftGenMap()
!!
!! DESCRIPTION
!!
!! This routine does not send the actual data, but uses communication
!! to generate the information about the communication pattern to be expected
!! when actual data movement happens.
!!
!!***

subroutine gr_pfftGenMap()
#include "Flash.h"
#include "constants.h"
#include "Pfft.h"

  use Logfile_interface, ONLY : Logfile_open, Logfile_close, Logfile_stamp
  use gr_pfftData, ONLY : pfft_procGrid, pfft_ndim, pfft_inLen, pfft_myPE, &
       pfft_me, pfft_comm, pfft_mode
  use gr_pfftReconfigData, ONLY : pfft_pencilSize, pfft_procLookup, &
       pfft_minRefLev, pfft_maxRefLev, pfft_oneRefLev, pfft_flashMsgList, &
       pfft_pencilMsgList, pfft_logMode, pfft_logUnit, &
       pfft_mpiTag1, pfft_mpiMetadataType
  use gr_pfftinterface, ONLY : gr_pfftGridPointTable, gr_pfftGenSingleMap
  use gr_pfftCommunicationPattern, ONLY : ut_mpiOverlapSendRecv, &
       common_state, send_state, recv_state
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_pfftReconfigFn, ONLY : gr_pfftFnArgSendMetadata, &
       gr_pfftFnArgRecvMetadata, gr_pfftFnArgDeliverMetadata
  use gr_pfftMessageList, ONLY : map_impure_fn

  implicit none
  include "Flash_mpi.h"
  interface
     integer function gr_pfft_allocate_buffers(item)
       use gr_pfftMessageNode, ONLY : message_node
       implicit none
       type(message_node), pointer :: item
     end function gr_pfft_allocate_buffers
  end interface

  type(common_state), pointer :: commonState
  type(send_state), pointer :: sendState
  type(recv_state), pointer :: recvState
  real :: secondsTillAbort
  integer, dimension(1:MDIM) :: gPencilStart, gPencilEnd
  integer :: commBufMemEstimate, err, axis, ierr
  logical :: sourceOnSingleLevelm, useBarriers
  character (len=*), parameter  :: logStr1d = "(a,1i5,a,1i5)"
  character (len=*), parameter  :: logStr2d = "(a,2i5,a,2i5)"
  character (len=*), parameter  :: logStr3d = "(a,3i5,a,3i5)"
#if NDIM == 1
    character (len=len(logStr1d)) :: logStr
    logStr = logStr1d
#elif NDIM == 2
    character (len=len(logStr2d)) :: logStr
    logStr = logStr2d
#elif NDIM == 3
    character (len=len(logStr3d)) :: logStr
    logStr = logStr3d
#endif

  if(pfft_ndim == 1) then
     return
  end if  

  !Stores values in pfft_procLookup.
  call gr_pfftGridPointTable(pfft_inLen)

  !Initialise any variables, arrays, datatypes.
  call gr_pfftInitMapData()
  

  !Calculate the number of grid points this processors will own in pencil space.
  !----------------------------------------------------------------------------
  do axis = 1, pfft_ndim
     gPencilStart(axis) = &
          pfft_procLookup(axis) % procInfo(pfft_me(axis)) % globalStartGridPoint
     gPencilEnd(axis) = &
          pfft_procLookup(axis) % procInfo(pfft_me(axis)) % globalEndGridPoint
  end do
  
  if ( (any(gPencilStart(1:pfft_ndim) == NONEXISTENT)) .or. &
       (any(gPencilEnd(1:pfft_ndim) == NONEXISTENT)) ) then
     !A NONEXISTENT means this processor will own no points in pencil space.
     !However, we still need this processor because it has FLASH grid points 
     !which must be sent to a pencil space processor.
     pfft_pencilSize = 0
     print *, "(INFO) Processor:", pfft_myPE, "has no pencil grid points."
  else
     if (pfft_logMode .eqv. .true.) then
        write(pfft_logUnit,logStr) " [gr_pfftGenMap]: "//&
             "My pencil grid points - start:", gPencilStart(1:NDIM), &
             "    end:", gPencilEnd(1:NDIM)
     end if
     pfft_pencilSize = product(&
          (gPencilEnd(1:pfft_ndim) - gPencilStart(1:pfft_ndim) + 1))
  end if
  !----------------------------------------------------------------------------


  !Generate a module-level map that is communicated using ut_mpiOverlapSendRecv.
  if (pfft_logMode .eqv. .true.) then
     write(pfft_logUnit,*) ""
     write(pfft_logUnit,*) "[gr_pfftGenMap]: Generating map 1"
     write(pfft_logUnit,*) "-------------------------------------------------------------"
  end if
  if (pfft_mode == PFFT_SINGLE_LEVEL) then
     call gr_pfftGenSingleMap(pfft_oneRefLev, leafMapMode=.false.)
  else if (pfft_mode == PFFT_MIXED_LEVEL) then
     call gr_pfftGenSingleMap(pfft_oneRefLev, leafMapMode=.true.)
  else
     call Driver_abortFlash("[gr_pfftGenMap]: Mode not recognised")
  end if


  if (pfft_logMode .eqv. .true.) then
     write(pfft_logUnit,*) ""
     write(pfft_logUnit,*) "[gr_pfftGenMap]: Communicating map 1"
     write(pfft_logUnit,*) "-------------------------------------------------------------"
  end if
  secondsTillAbort = 30.0
  useBarriers = .true.

  allocate(commonState, sendState, recvState, STAT=err)
  if (err /= 0) then
     call Driver_abortFlash("[gr_pfftGenMap]: State allocation failed.")
  end if

  commonState % mpiTag = pfft_mpiTag1 
  commonState % mpiType = pfft_mpiMetadataType
  commonState % mpiComm = pfft_comm(IAXIS)
  commonState % direction = TO_PFFT

  sendState % msgList => pfft_flashMsgList
  nullify (sendState % activeMsg)
  sendState % numMsgSent = 0

  recvState % msgList => pfft_pencilMsgList
  nullify (recvState % activeMsg)
  recvState % numMsgRecv = 0
  recvState % numGridPointsRecv = 0
  recvState % totGridPoints = pfft_pencilSize

  call ut_mpiOverlapSendRecv &
       (gr_pfftFnArgSendMetadata, gr_pfftFnArgRecvMetadata, &
       gr_pfftFnArgDeliverMetadata, commonState, sendState, recvState, &
       secondsTillAbort, useBarriers)

  !Double check whether I should nullify sendState % msgList
  nullify (sendState % activeMsg, sendState % msgList, recvState % msgList)
  deallocate(commonState, sendState, recvState, STAT=err)
  if (err /= 0) then
     call Driver_abortFlash("[gr_pfftGenMap]: State deallocation failed.")
  end if


  
  !At this point all processors have the LEAF metadata they need.
  !But we need another round of communication to discover the 
  !processors that actually own the FLASH blocks at level=solveLevel.
  !This is because the last communication round included restricted 
  !LEAF block data (where restriction done in place).  However, the 
  !FLASH block covering the same space at the restricted level may 
  !exist on another processor (needed for return journey PFFT -> FLASH).
  if (pfft_mode == PFFT_MIXED_LEVEL) then
     !Execute when the LEAF blocks are at different refinement levels.
     if ( (pfft_minRefLev /= pfft_maxRefLev) .or. &
          ((pfft_minRefLev == pfft_maxRefLev) .and.&
          (pfft_oneRefLev /= pfft_maxRefLev)) ) then

        call Driver_abortFlash("[gr_pfftGenMap]: Not yet coded!")

        !We need to generate a map for the blocks at "solveLevel".
        call gr_pfftGenSingleMap(pfft_oneRefLev)

        !call ut_mpiOverlapSendRecv()
        
        !The processors that sent non-LEAF solveLevel block fragments 
        !must now post a receive for each fragment.  In the sent message 
        !we inform the coarse level block about what data they will 
        !receive and from what fragment.
        !call gr_pfftNotifySolveLevelBlock()

     end if
  end if


  !Each processor has constructed a list of FLASH nodes and a list 
  !of PFFT nodes.  Now we traverse the list and allocate the data buffer
  !within each message node.
  call map_impure_fn(gr_pfft_allocate_buffers, pfft_flashMsgList)
  call map_impure_fn(gr_pfft_allocate_buffers, pfft_pencilMsgList)
         
end subroutine gr_pfftGenMap


integer function gr_pfft_allocate_buffers(item)
  use gr_pfftMessageNode, ONLY : message_node
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  type(message_node), pointer :: item
  integer :: err
  allocate (item % buf (item % bufSize), STAT=err)
  if (err /= 0) then
     call Driver_abortFlash("[gr_pfft_allocate_buffers:" //&
          "Allocation failed")
  end if
  gr_pfft_allocate_buffers = 0
end function gr_pfft_allocate_buffers
