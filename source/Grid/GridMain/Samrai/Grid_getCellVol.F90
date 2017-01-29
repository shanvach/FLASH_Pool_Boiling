!!****if* source/Grid/GridMain/Samrai/Grid_getCellVol
!!
!! NAME
!!  Grid_getCellVol
!!
!! SYNOPSIS
!!  #include "Flash.h"
!!  #include "constants.h"
!!
!!  Grid_getCellVol(integer(IN) :: axis(MDIM),
!!                  integer(IN) :: blockid,
!!                  logical(IN) :: guardcell, 
!!                  real(OUT)   :: cellVolume(size(1), size(2), size(3)),
!!                  integer(IN) :: size(MDIM))
!!  
!! DESCRIPTION 
!!  
!!  Gets cell volumes for the specified cells in a given block.  This routine
!!  can return a single cell's volume, a row of cells' volumes, a plane of cells'
!!  volumes or the cell volumes for the entire block. 
!!  
!! ARGUMENTS 
!!
!!  axis(MDIM) specifies the integer index coordinates of the cells
!!             A WILDCARD along any of the axes tells the funtion
!!             to get all the cell volumes along that dimension. 
!!             If the axis doesn't exist (for example KAXIS in a 2D 
!!             problem, its value defaults to KHI_GC, which should be 1).
!!
!!             WILDCARD is a defined constant located in "constants.h"  
!!             
!!             IF zero axes is WILDCARD, a single cell volume is returned
!!             IF only one axis is WILDCARD, a row of cell volumes is returned
!!             IF any two axes are WILDCARD, a plane of cell volumes is returned
!!             IF all three axes = WILDCARD, all cell volumes in block are returned
!!
!!             ex. to put only one cell:
!!             if axis(1) = i, axis(2) = j and axis(3) = k then the routine
!!             will get the cell volumes for a the particular cell(i,j,k) in a given
!!             block
!!
!!             ex. to put a row:
!!             if axis(1) = WILDCARD, axis(2) = j and axis(3) = k
!!             all cell volumes in the ith direction, with position j and k
!!             will be returned
!!
!!             ex. to put an i,j plane: 
!!             if axis(1) = WILDCARD, axis(2) = WILDCARD and axis(3) = k
!!             all cell vollumes in the ith and jth direction, at position k
!!             will be returned
!! 
!!  blockid - integer local blockid
!!
!!  guardcell - logical value, if true, all cells, interior and guardcells are returned.  
!!              If false only the interior cells are returned
!!              
!!
!!  cellVolume - 3 dimensional real array returned holding the cell volumes.  
!!               Each dimension of cellVolume holds the cell volumes of one axis of
!!               a block.  The particular dimensions are depend on if the argument 
!!               guardcell is true or false. The compiler needs to know the size of
!!               cellVolume so we explicitly send in the sizes in the next argument
!!               "size" see below
!!
!!  size(MDIM) -  integer array specifying the dimensions for argument cellVolume
!!          
!!          The values of size depend of if guardcells are included or not
!!          size(1) holds the number of cells in the i direction
!!          IF (guardcell==.true.) then size(1)= GRID_IHI_GC, ELSE size(1)= NXB
!!
!!          size(2) holds the number of cells in the j direction
!!          IF (guardcell==.true.) then size(2)= GRID_JHI_GC, ELSE size(2)= NXB
!!
!!          size(3) holds the number of cells in the k direction
!!          IF (guardcell==.true.) then size(3)= GRID_KHI_GC, ELSE size(3)= NXB
!!
!!
!!
!! EXAMPLE
!!
!!  usage of axis argument
!!                 
!!                 <WILDCARD,2> in axis
!!                 for a 2D problem tell the subroutine to fetch volumes of cells
!!                 in the second column.
!!                 Similarly <WILDCARD,WILDCARD,4> in 3D fetches the 4th XY plane cell volumes.
!!                 If all the values in axis are WILDCARDS, then the entire block is fetched.
!!
!!
!! NOTES
!!  if the constant WILDCARD is passed in as an argument the calling routine needs
!!  to include the header file "constants.h"
!! 
!!  IAXIS, JAXIS and KAXIS (defined as 1,2,and3) are used in this routine to increase 
!!  readability of the code.  They are also defined in "constants.h" 
!! 
!!
!!***

subroutine Grid_getCellVol(axis, blockid, guardcell, cellvolume ,size)

  use Grid_data, ONLY : gr_maxPatches

  implicit none

#include "constants.h"

  integer, dimension(MDIM), intent(in) :: axis, size
  integer, intent(in) :: blockid
  logical, intent(in) :: guardcell
  real, dimension(size(1), size(2), size(3)), intent(out) :: cellvolume

  integer :: levelNum, patchNum
  
  !translate blockID to patchNum and levelNum
  levelNum = blockID / gr_maxPatches

  patchNum = mod(blockID, gr_maxPatches)

  
  !verify if logicals can get translated correctly to booleans in C++
  call samrai_get_cell_volume(axis, patchNum, levelNum, guardcell)

  return
end subroutine Grid_getCellVol
