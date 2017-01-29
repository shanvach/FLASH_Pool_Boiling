!!****if* source/Grid/GridMain/Samrai/Grid_getLocalNumBlks
!!
!! NAME
!!  Grid_getLocalNumBlks
!!
!! SYNOPSIS
!!  #include "Flash.h"
!!
!!  Grid_getLocalNumBlks(integer(OUT) :: numBlocks)
!!  
!! DESCRIPTION 
!!  Get the number of local blocks on a processor 
!!
!!***


subroutine Grid_getLocalNumBlks(numBlocks)

implicit none
  integer,intent(out) :: numBlocks

  call samrai_get_local_num_patches(numBlocks)

  return
end subroutine Grid_getLocalNumBlks
