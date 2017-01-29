!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/ChooseLevel/gr_pfftFragmentNode
!!
!! NAME
!!  gr_pfftFragmentNode
!!
!! SYNOPSIS
!!
!!  use gr_pfftFragmentNode
!!
!! DESCRIPTION
!!
!!  Object which represents a block fragment.
!!   
!!***
module gr_pfftFragmentNode
#include "constants.h"
#include "Flash.h"
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  include "Flash_mpi.h"

  integer, save :: m_unitNum, m_countID
  logical, save :: m_logEvents

  type node_metadata
     integer :: srcFlashProcID      !Process containing source data.
     integer :: srcFlashBlockID     !Block containing source data.
     integer :: srcFlashBlockType   !Node type of source block.
     integer :: srcFlashBlockRefLev !Refinement level of source block.
     integer :: solnFlashProcID     !Process containing solution data.
     integer :: solnFlashBlockID    !Block containing solution data.
     integer :: pfftProcID          !PFFT process.
     integer :: bufStart            !Data start position in message node.
     integer :: bufSize             !Data size in message node.
     
     !Block fragment coordinates (solve level).
     integer, dimension(1:MDIM) :: srcFlashStartPos, srcFlashEndPos

     !Block fragment coordinates (source block level).
     integer, dimension(1:MDIM) :: srcFlashActualStartPos, &
          srcFlashActualEndPos
     
     !Block fragment coordinates in the PFFT pencil.
     integer, dimension(1:MDIM) :: pfftStartPos, pfftEndPos
  end type node_metadata


  type fragment_node
     type(node_metadata) :: metadata
     type(fragment_node), pointer :: next, prev
     integer, dimension(MPI_STATUS_SIZE) :: metadata_status
     integer :: metadata_request, ID
  end type fragment_node

contains 

  !The user must provide: create_node, destroy_node and print_node.
  subroutine create_node(item)
    implicit none
    type(fragment_node), pointer :: item
    integer :: err

    allocate(item, STAT=err)
    if (err /= 0) then
       call Driver_abortFlash ("[gr_pfftFragmentNode]: "//&
            "Error in node allocation")
    end if    
    nullify(item % next, item % prev)

    item % ID = m_countID
    item % metadata % bufStart = 1
    item % metadata % bufSize = 0
    if (m_logEvents .eqv. .true.) then
       write(m_unitNum,*) "[gr_pfftFragmentNode]: "//&
            "Allocated fragment node ", item % ID
    end if
    m_countID = m_countID + 1
  end subroutine create_node


  subroutine destroy_node(item)
    implicit none
    type(fragment_node), pointer  :: item
    integer :: err, ID

    ID = item % ID
    deallocate(item, STAT=err)
    if (err /= 0) then
       call Driver_abortFlash ("[gr_pfftFragmentNode]: "//&
            "Error in node deallocation")
    end if
    nullify(item)

    if (m_logEvents .eqv. .true.) then
       write(m_unitNum,*) "[gr_pfftFragmentNode]: "//&
            "Freed fragment node ", ID
    end if
  end subroutine destroy_node


  subroutine print_node(unitNumber, item)    
    implicit none
    integer, intent(IN) :: unitNumber
    type(fragment_node), pointer :: item

    character (len=*), parameter  :: baseStr = &
         "(/ a,i5 / a,i5 / a,i5 / a,i5 / a,i5 / a,i5 / a,i5 / a,i5 / a,i10 / a,i10 /"
    character (len=*), parameter  :: logStr1d = &
         "a,1i5 / a,1i5 / a,1i5 / a,1i5 / a,1i5 / a,1i5 /)"
    character (len=*), parameter  :: logStr2d = &
         "a,2i5 / a,2i5 / a,2i5 / a,2i5 / a,2i5 / a,2i5 /)"
    character (len=*), parameter  :: logStr3d = &
         "a,3i5 / a,3i5 / a,3i5 / a,3i5 / a,3i5 / a,3i5 /)"
#if NDIM == 1
    character (len=len(baseStr)+len(logStr1d)) :: logStr
    logStr = baseStr // logStr1d
#elif NDIM == 2
    character (len=len(baseStr)+len(logStr2d)) :: logStr
    logStr = baseStr // logStr2d
#elif NDIM == 3
    character (len=len(baseStr)+len(logStr3d)) :: logStr
    logStr = baseStr // logStr3d
#endif

    write(unitNumber,logStr) &
         " [gr_pfftFragmentNode]: Fragment node local ID:", item % ID, &
         "    srcFlashProcID:", item % metadata % srcFlashProcID, &
         "    srcFlashBlockID:", item % metadata % srcFlashBlockID, & 
         "    srcFlashBlockType:", item % metadata % srcFlashBlockType, &
         "    srcFlashBlockRefLev:", item % metadata % srcFlashBlockRefLev, &
         "    solnFlashProcID:", item % metadata % solnFlashProcID, &
         "    solnFlashBlockID:", item % metadata % solnFlashBlockID, &
         "    pfftProcID:", item % metadata % pfftProcID, &
         "    bufStart:", item % metadata % bufStart, &
         "    bufSize:", item % metadata % bufSize, &
         "    srcFlashStartPos:", item % metadata % srcFlashStartPos(1:NDIM), &
         "    srcFlashEndPos:", item % metadata % srcFlashEndPos(1:NDIM), &
         "    srcFlashActualStartPos:", item % metadata % srcFlashActualStartPos(1:NDIM), &
         "    srcFlashActualEndPos:", item % metadata % srcFlashActualEndPos(1:NDIM), &
         "    pfftStartPos:", item % metadata % pfftStartPos(1:NDIM), &
         "    pfftEndPos:", item % metadata % pfftEndPos(1:NDIM)
    !call flush(unitNumber)
  end subroutine print_node

end module gr_pfftFragmentNode
