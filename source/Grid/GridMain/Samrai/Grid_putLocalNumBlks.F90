!!****if* source/Grid/GridMain/Samrai/Grid_putLocalNumBlks
!!
!! NAME
!!  Grid_putLocalNumBlks
!!
!! SYNOPSIS
!!
!!  Grid_putLocalNumBlks(integer(IN) :: numBlocks)
!!  
!! DESCRIPTION 
!!  Put the number of local blocks on a processor 
!!  Only used in restart capabilities
!!
!!***

subroutine Grid_putLocalNumBlks(numBlocks)
implicit none
  integer,intent(in) :: numBlocks

  return
end subroutine Grid_putLocalNumBlks
