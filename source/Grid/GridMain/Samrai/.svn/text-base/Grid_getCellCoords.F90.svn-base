!!****if* source/Grid/GridMain/Samrai/Grid_getCellCoords
!!
!! NAME
!!  Grid_getCellCoords
!!
!! SYNOPSIS
!!  #include "Flash.h"
!!
!!  Grid_getCellCoords(integer(IN)  :: axis,
!!                      integer(IN) :: blockid, 
!!                      integer(IN) :: edge, 
!!                      logical(IN) :: guardcell, 
!!                      real(OUT)   :: coordinates(size),
!!                      integer(IN) :: size)
!!  
!! DESCRIPTION
!!    This subroutine is an accessor function that gets the coordinates of
!!    the cells in a given block.
!!    Coordinates are retrieved one axis at a time, meaning you can get the i, j, _or_ k
!!    coordinates with one call.  If you want all the coordinates, all axes, you
!!    need to call Grid_getCellCoords 3 times, one for each axis.
!!
!!
!!
!! ARGUMENTS
!!            
!!   axis - specifies the integer index coordinates of the cells being retrieved.
!!          axis can have one of three different values, IAXIS, JAXIS or KAXIS 
!!          (defined in constants.h as 1,2 and 3)
!!
!!   blockId - integer block number
!!
!!   edge - integer value with one of four values, LEFT, RIGHT, CENTER or WILDCARD
!!          (These are #defined constants located in constants.h)
!!        The edge argument specifies what side of the zone to get, the CENTER
!!        point, the LEFT side or the RIGHT side of the zone.  A user can get all
!!        three edges, LEFT, CENTER and RIGHT coordinates of a zone with the
!!        WILDCARD flag.
!!
!!
!!   guardcell - logical value, if true, coordinates of a block including guardcells
!!                are returned, if false, only interior cell coords are returned
!!
!!          
!!   coordinates - The array holding the data returning the coordinate values
!!             coordinates is of size (see below)
!!           
!!   size : a 2d integer array specifying the dimensions for coordinates
!!        
!!          The values of the 'size' array depend on other arguments.
!!   
!!         size = size of the vector returned.  If guardcell true (and axis =IAXIS) 
!!               then
!!                   size(2) = nguard*2 + nxb for a fixed block size example.
!!                   if guardcell = false then size(2) = nxb
!!
!!
!!
!!***


subroutine Grid_getCellCoords(axis, blockid, edge, guardcell, coordinates, size)
  
  use Grid_data, ONLY : gr_maxPatches

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(in) :: axis, blockid, edge
  integer, intent(in) :: size
  logical, intent(in) :: guardcell
  real,intent(out), dimension(size) :: coordinates
  

  integer :: levelNum, patchNum, i
  real, dimension(MDIM) :: singleCoord, del

  !translate blockID to patchNum and levelNum
  levelNum = blockID / gr_maxPatches

  patchNum = mod(blockID, gr_maxPatches)

  !del returns the grid spacing
  call samrai_get_cell_coords(patchNum, levelNum, guardcell, singleCoord, del)

  
  !we are going to calculate the coordinates ourselves
  !katie using, del return array of coordinates
  do i=0, size -1
     coordinates(i) = singleCoord(axis) + i*del(axis) + del(axis)/2
  end do


end subroutine Grid_getCellCoords
