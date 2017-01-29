!!****if* source/Grid/GridMain/Samrai/Grid_getBlkBC
!!
!! NAME
!!  Grid_getBlkBC
!!
!! SYNOPSIS
!!  #include "Flash.h"
!!
!!  Grid_getBlkBC(integer(IN)  :: blockId,
!!                        integer(OUT) :: faces(2,MDIM))
!!                    
!! DESCRIPTION 
!!  For each face of the block, this function finds out if it is
!!  on the physical boundary, and if so, returns the boundary condition
!!  Returns the boundary conditions. 
!!   
!! ARGUMENTS 
!!
!!  blockId - the local blockId, 
!!  faces   - the array to return the face related information in
!!
!!            the first index of the array can take on values LOW or
!!            HIGH, and the second index can be IAXIS, JAXIS or KAXIS
!!            The returned value for any entry is the array is NO_BOUDARY
!!            if that face along that axis is not on boundary, otherise
!!            it contains the boundary condition
!!
!!***


subroutine Grid_getBlkBC(blockId, faces)

  use Grid_data, ONLY : gr_maxPatches

implicit none
#include "constants.h"

  integer, intent(in) :: blockId
  integer, dimension(2,MDIM),intent(out):: faces

  
  integer :: levelNum, patchNum
  !translate blockID to patchNum and levelNum

  levelNum = blockID / gr_maxPatches

  patchNum = mod(blockID, gr_maxPatches)

  call samrai_get_patch_bc(patchNum, levelNum, faces)

  return
end subroutine Grid_getBlkBC
