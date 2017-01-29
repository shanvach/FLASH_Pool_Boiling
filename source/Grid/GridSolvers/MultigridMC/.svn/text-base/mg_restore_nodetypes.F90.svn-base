!*******************************************************************************

!  Routine:     mg_restore_nodetypes

!  Description:  restores the nodetypes and newchild arrays to pre-multigrid states    
!                calls the expensive get_tree_nodetypes
!               



subroutine mg_restore_nodetypes (level)

!===============================================================================


  use mg_common, ONLY: nodetype_save, newchild_save
 
  use Grid_data, ONLY : gr_meshMe,gr_meshNumProcs

  use tree, only : nodetype,newchild

  use Grid_interface,    ONLY : Grid_getLocalNumBlks
 
  implicit none

  integer, intent(in) :: level

  integer :: lnblocks, lb

  
  call Grid_getLocalNumBlks(lnblocks)

  do lb = 1, lnblocks
     nodetype(lb) = nodetype_save(lb)
     newchild(lb) = newchild_save(lb)
  enddo

!!$  call mpi_amr_read_guard_comm_mg(numPEs, level)
!!$  call mpi_amr_read_prol_comm_mg(numPEs,  level)
!!$  call mpi_amr_read_flux_comm_mg(numPEs,  level)
!!$  call mpi_amr_read_restrict_comm_mg(numPEs,  level)

  call amr_get_new_nodetypes (gr_meshNumProcs, gr_meshMe, level)

end subroutine mg_restore_nodetypes

