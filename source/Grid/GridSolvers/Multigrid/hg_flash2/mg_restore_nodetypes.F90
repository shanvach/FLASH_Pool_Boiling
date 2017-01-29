!*******************************************************************************

!  Routine:     mg_restore_nodetypes

!  Description:  restores the nodetypes and newchild arrays to pre-multigrid states    
!                calls the expensive get_tree_nodetypes
!               



subroutine mg_restore_nodetypes ()

!===============================================================================

  use perfmon
  use mg_common, ONLY: nodetype_save, newchild_save
  use dBase, ONLY: dBaseTreePtrNodeType, &
       dBasePropertyInteger, &
       dBaseTreePtrNewChild
  
  implicit none

  integer, save :: myPE, numPEs 
  integer :: lnblocks, lb
  integer, pointer, dimension(:), save :: nodetype
  logical, pointer, dimension(:), save :: newchild
  logical, save :: first_call = .true.

!===============================================================================
  
  if (first_call) then
     first_call = .false.
     myPE = dBasePropertyInteger("MyProcessor")
     numPEs = dBasePropertyInteger("NumberOfProcessors")
     nodetype =>dBaseTreePtrNodeType()
     newchild =>dBaseTreePtrNewChild()
  end if
  
  lnblocks = dBasePropertyInteger("LocalNumberOfBlocks")

  do lb = 1, lnblocks
     nodetype(lb) = nodetype_save(lb)
     newchild(lb) = newchild_save(lb)
  enddo
!  call timer_start("get_tree_nodetypes")
  call get_tree_nodetypes (NumPEs, MyPE)
!  call timer_stop("get_tree_nodetypes")

!===============================================================================

return
end subroutine mg_restore_nodetypes

