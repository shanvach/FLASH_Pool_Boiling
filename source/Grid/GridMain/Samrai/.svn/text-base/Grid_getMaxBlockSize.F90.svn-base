!!****if* source/Grid/GridMain/Samrai/Grid_getMaxBlockSize
!!
!! NAME
!!  Grid_getMaxBlockSize
!!
!! SYNOPSIS
!!
!!  Grid_getMaxBlockSize(integer(OUT) :: maxBlockSize(MDIM))
!!  
!! DESCRIPTION 
!!  In non fixedblocksize mode, get the maximum allowed block size to 
!!  allocate a single block array
!!
!!***

subroutine Grid_getMaxBlockSize(maxBlockSize)
implicit none
# include "constants.h"  

  integer,dimension(MDIM), intent(out) :: maxBlockSize

  call samrai_get_max_patch_size(maxBlockSize)

end subroutine Grid_getMaxBlockSize

