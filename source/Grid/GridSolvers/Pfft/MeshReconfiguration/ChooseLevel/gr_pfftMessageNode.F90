!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/ChooseLevel/gr_pfftMessageNode
!!
!! NAME
!!  gr_pfftMessageNode
!!
!! SYNOPSIS
!!
!!  use gr_pfftMessageNode
!!
!! DESCRIPTION
!!
!!  Derived data type which represents a message.
!!   
!!***
module gr_pfftMessageNode
#include "constants.h"
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_pfftFragmentList, ONLY : fragment_list, &
       initialise_fragment_list => initialise_list, &
       finalise_fragment_list => finalise_list
  implicit none
  include "Flash_mpi.h"

  !We can't refer to pfft_unitNum, pfft_logMode because these 
  !are in a module that uses a node datatype from this module.
  !This would create a circular-dependency. 
  integer, save :: m_unitNum, m_countID
  logical, save :: m_logEvents

  type message_node
     !A message node contains a list describing the various block fragments.
     type(fragment_list), pointer :: fragmentList
     type(message_node), pointer :: next, prev

     !Buffer containing grid data (Must be a pointer as in a type).
     real, pointer, dimension(:) :: buf

     integer, dimension(MPI_STATUS_SIZE) :: data_status
     integer :: bufSize, data_request, PE_partner, ID
     logical :: msgDelivered
  end type message_node

contains

  !The user must provide: create_node, destroy_node and print_node.
  subroutine create_node(item)
    implicit none
    type(message_node), pointer :: item
    integer :: err

    allocate(item, STAT=err)
    if (err /= 0) then
       call Driver_abortFlash ("[gr_pfftMessageNode]: "//&
            "Error in node allocation")
    end if
    nullify(item % next, item % prev, item % buf, item % fragmentList)

    item % ID = m_countID
    item % bufSize = 0
    item % PE_Partner = -1
    item % msgDelivered = .false.
    call initialise_fragment_list(item % fragmentList)
    if (m_logEvents .eqv. .true.) then
       write(m_unitNum,*) "[gr_pfftMessageNode]: "//&
            "Allocated message node ", item % ID
    end if
    m_countID = m_countID + 1
  end subroutine create_node


  !We have a custom deallocation subroutine because nodes may 
  !contain lots of pieces of data which must be deallocated appropriately.
  subroutine destroy_node(item)
    implicit none
    type(message_node), pointer  :: item
    integer :: err, ID

    ID = item % ID
    call finalise_fragment_list(item % fragmentList)
    nullify(item % fragmentList)

    if (associated(item % buf)) then
       deallocate(item % buf, STAT=err)
       if (err /= 0) then
          call Driver_abortFlash ("[gr_pfftMessageNode]: "//&
               "Error in node % buf deallocation")
       end if
       nullify(item % buf)
    end if

    nullify(item % next, item % prev)
    deallocate(item, STAT=err)
    if (err /= 0) then
       call Driver_abortFlash ("[gr_pfftMessageNode]: "//&
            "Error in node deallocation")
    end if
    nullify(item)

    if (m_logEvents .eqv. .true.) then
       write(m_unitNum,*) "[gr_pfftMessageNode]: "//&
            "Freed message node ", ID
    end if
  end subroutine destroy_node


  subroutine print_node(unitNumber, item)    
    implicit none
    integer, intent(IN) :: unitNumber
    type(message_node), pointer :: item
    character (len=*), parameter  :: logStr = &
         "(/ a,i4 / a,i4 / a,i10 /)"
    write(unitNumber,logStr) &
         " [gr_pfftMessageNode]: Message node local ID", item % ID, &
         "    PE_partner:", item % PE_partner, & 
         "    bufSize:", item % bufSize
    !call flush(unitNumber)
  end subroutine print_node

end module gr_pfftMessageNode
