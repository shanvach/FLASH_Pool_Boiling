subroutine gr_pfftGetMessageNode(PE, messageList, messageNode)
  use gr_pfftMessageNode, ONLY : message_node, create_node
  use gr_pfftMessageList, ONLY : message_list, push_back, find_matching_node
  use gr_pfftReconfigData, ONLY : pfft_logMode, pfft_logUnit
  implicit none

  integer, intent(IN) :: PE
  type(message_list), pointer :: messageList !May be extended in this subroutine.
  type(message_node), pointer :: messageNode !We always return a valid messageNode.

  interface 
     logical function gr_pfftFnArgMatchPE(messageNode, PE)
       use gr_pfftMessageNode, ONLY : message_node
       implicit none
       type(message_node), pointer :: messageNode
       integer, intent(IN) :: PE
     end function gr_pfftFnArgMatchPE
  end interface
  nullify(messageNode)

  !The first thing we need to do is identify whether we already have 
  !a message node assigned to PE.  If we do, get a pointer 
  !to this node, otherwise create a new node.  
  call find_matching_node(gr_pfftFnArgMatchPE, PE, &
       messageList, messageNode)

  if (.not.associated(messageNode)) then
     call create_node(messageNode)
     messageNode % PE_partner = PE
     call push_back(messageList, messageNode)
     if (pfft_logMode .eqv. .true.) then
        write(pfft_logUnit,*) "[gr_pfftGetMessageNode]: "//&
             "Created new message node for process ", PE
     end if
  else
     if (pfft_logMode .eqv. .true.) then
        write(pfft_logUnit,*) "[gr_pfftGetMessageNode]: "//&
             "Retrieved message node for process ", PE
     end if
  end if
end subroutine gr_pfftGetMessageNode


logical function gr_pfftFnArgMatchPE(messageNode, PE)
  use gr_pfftMessageNode, ONLY : message_node
  implicit none
  type(message_node), pointer :: messageNode
  integer, intent(IN) :: PE

  if (messageNode % PE_partner == PE) then
     gr_pfftFnArgMatchPE = .true.
  else
     gr_pfftFnArgMatchPE = .false.
  end if
end function gr_pfftFnArgMatchPE
