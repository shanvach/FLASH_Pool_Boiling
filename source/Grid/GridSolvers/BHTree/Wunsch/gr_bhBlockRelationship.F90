!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhBlockRelationship
!!
!! NAME
!!
!!  gr_bhBlockRelationship
!!
!!
!! SYNOPSIS
!!
!!  gr_bhBlockRelationship(
!!        integer, intent(in) :: block,
!!        integer, intent(in) :: tr,
!!        integer, intent(in) :: cpu
!!        )
!!
!! DESCRIPTION
!!
!!  Determines if two blocks are neighbours (or if it's the same block) and 
!!  the direction in which they neighbour
!!
!! ARGUMENTS
!!
!!  block   - ID of a block for which neighbours are tested
!!  tr      - ID of a block whose potential neighbourship is tested
!!  cpu     - cpu of a block whose potential neighbourship is tested
!!
!!***

integer function gr_bhBlockRelationship(block, tr, cpu)
  use gr_bhData, ONLY : gr_bhTreeSurbox, gr_bhTreeMyPE

  implicit none
  integer, intent(in) :: block, tr, cpu
  integer :: temp

  gr_bhBlockRelationship = -1
  do temp = 1,27
    if ( (tr .eq. gr_bhTreeSurbox(1, temp, block, gr_bhTreeMyPE)) .and. &
         (cpu .eq. gr_bhTreeSurbox(2, temp, block, gr_bhTreeMyPE))) then
      gr_bhBlockRelationship = temp
    endif
  enddo

  return
end function gr_bhBlockRelationship


