!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhGetTreeSize
!!
!! NAME
!!
!!  gr_bhGetTreeSize
!!
!!
!! SYNOPSIS
!!
!!   gr_bhGetTreeSize(
!!        integer,intent(in) :: level
!!        )
!!
!! DESCRIPTION
!!
!!   Determines the size of the tree of a given level.
!!
!! ARGUMENTS
!!
!!   level - level of the tree of which the size is calculated
!!
!!***



integer function gr_bhGetTreeSize(level)

  use gr_bhData, ONLY : gr_bhTreeLevels, gr_bhTreeBS
  implicit none
  integer,intent(in) :: level

  if (level == gr_bhTreeLevels) then
    gr_bhGetTreeSize = 8**gr_bhTreeLevels + 4*(8**gr_bhTreeLevels - 1)/7
  else
    gr_bhGetTreeSize = 4*(8**(level+1) - 1)/7
  endif

  return
end function gr_bhGetTreeSize

