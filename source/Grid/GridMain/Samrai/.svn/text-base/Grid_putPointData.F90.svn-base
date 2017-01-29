!!****if* source/Grid/GridMain/Samrai/Grid_putPointData
!!
!! NAME
!!  Grid_putPointData
!!
!! SYNOPSIS
!!
!!  Grid_putPointData(integer(IN) :: blockID,
!!                    integer(IN) :: cellType,
!!                    integer(IN) :: variable,
!!                    integer(IN) :: beginCount, 
!!                    integer(IN) :: position(MDIM),
!!                    real(IN)    :: datablock)
!!                    
!!  
!! DESCRIPTION 
!!  
!!  Puts a point of simulation data for specified variable and cells.
!!  
!!  The user is
!!  allowed to specify if indice counting should begin at the exterior edge
!!  of the block, (that is including guardcells)
!!  or the interior edge of the block 
!!
!!  
!! ARGUMENTS 
!!
!!  blockid : the local blockid
!!
!!  cellType : integer value specifying the type of cell data.  The options are
!!             CENTER, FACEX, FACEY, FACEZ defined in constants.h indicated
!!             whether the data is cell centered of face data. 
!!
!!  variable : integer value that specifies which variable to put into storage.
!!             ex: DENS_VAR, PRES_VAR as defined in Flash.h 
!!
!!  
!!
!!  beginCount : tells the routine where to start index counting.  beginCount can
!!               be set to INTERIOR or EXTERIOR.  If INTERIOR is specified
!!               guardcell indices are not included and index 1 is the first interior cell. 
!!               If EXTERIOR is specified
!!               the first index, 1, is the left most guardcell.  See example
!!               below for more explanation.  (For most of the FLASH architecture code,
!!               we use EXTERIOR.  Some physics routines, however, find it helpful 
!!               only to work on the internal parts of the blocks (without
!!               guardcells) and wish to keep loop indicies  
!!               going from 1 to NXB without having to worry about finding 
!!               the correct offset for the number of guardcells.) 
!!               (INTERIOR and EXTERIOR are defined in constants.h)
!!
!!
!!  position(MDIM):
!!           specifies the point to return
!!   
!!           position(1) = i
!!           position(2) = j
!!           position(3) = k
!!
!!           If a problem is only 2d, position(3) is ignored.  For 1d problems
!!           position(2) and position(3) are ignored.
!!
!!
!!  datablock : a real value containing the data
!!
!!
!! EXAMPLE  
!!  
!!    Here is a 3d block example putting a point of data, 
!!    for each block on a local processor.  We will put the point at 
!!    position i=5, j=6 and k=7.  beginCount is set to EXTERIOR
!!    meaning that the first cell of the block including guardcells is 
!!    index = 1.  If in this example, NGUARD = 4 then position 5 is the
!!    first interior cell in a dimension.
!!
!! 
!!    (very hard to draw this, especially the j axis. 
!!     picture is just meant to help show where 
!!     counting begins when "beginCount" is set to EXTERIOR)
!!    
!!     j
!!    
!!    16         - - - - - - - -  
!!    15         - - - - - - - - 
!!    14         - - - - - - - - 
!!    13         - - - - - - - - 
!!    12 - - - -|-|-|-|-|-|-|-|-|- - - -
!!    11 - - - -|-|-|-|-|-|-|-|-|- - - -
!!    10 - - - -|-|-|-|-|-|-|-|-|- - - -
!!     9 - - - -|-|-|-|-|-|-|-|-|- - - -
!!     8 - - - -|-|-|-|-|-|-|-|-|- - - -
!!     7 - - - -|-|-|-|-|-|-|-|-|- - - -
!!     6 - - - -|-|-|-|-|-|-|-|-|- - - -
!!     5 - - - -|-|-|-|-|-|-|-|-|- - - -
!!     4         - - - - - - - - 
!!     3         - - - - - - - - 
!!     2         - - - - - - - - 
!!     1         - - - - - - - - 
!! i     1 2 3 4 5 6 7 8 9 10111213141516 !!
!!
!!
!!
!!    #include "Flash.h"
!!    #include "constants.h"
!!
!!    ...
!!       
!!      integer ::    position(MDIM)
!!      integer ::    blockID
!!      real    ::    dataBlock
!!       
!!          position(1) = 5    
!!          position(2) = 6
!!          position(3) = 7
!!
!!
!!
!!          do blockID = 1, localNumBlocks
!!  
!!             call Grid_putPointData(blockID, CENTER, DENS_VAR, EXTERIOR, &
!!                               position, dataBlock)
!!  
!!          end do
!!
!!
!!  
!!    In this 2d block example we will put a point of data for the pressure variable
!!    for each block on a local processor.
!!    beginCount is set to INTERIOR, meaning that all the position indices
!!    will start where index 1 is the first interior cell of the block.
!!    In this example we will put a point where i=4, j=5
!!
!!    (hard to draw, but this is the idea, stars (*) are the cells to return
!!     notice the where indice counting starts when beginCount is set to INTERIOR)
!!            - - - - - - - - 
!!            - - - - - - - - 
!!            - - - - - - - - 
!!     j      - - - - - - - - 
!!     8 ----|-|-|-|-|-|-|-|-|----
!!     7 ----|-|-|-|-|-|-|-|-|----
!!     6 ----|-|-|-|-|-|-|-|-|----
!!     5 ----|-|-|-|-|-|-|-|-|----
!!     4 ----|-|-|-|*|-|-|-|-|----
!!     3 ----|-|-|-|-|-|-|-|-|----
!!     2 ----|-|-|-|-|-|-|-|-|----
!!     1 ----|-|-|-|-|-|-|-|-|----
!!            - - - - - - - - 
!!            - - - - - - - - 
!!            - - - - - - - - 
!!            - - - - - - - - 
!!         i  1-2-3-4 5-6-7-8 
!!
!!
!! 
!!    #include "Flash.h"
!!    #include "constants.h"
!!
!!    ...
!!       
!!      integer ::    position(MDIM)
!!      integer ::    blockID
!!      real    ::    dataBlock
!!       
!!          position(1) = 4    
!!          position(2) = 5
!!          position(3) = 1 !ignored since only 2d problem
!!
!!
!!
!!          do blockID = 1, localNumBlocks
!!  
!!             call Grid_putPointData(blockID, CENTER, PRES_VAR, INTERIOR, &
!!                               position, dataBlock)
!!  
!!          end do
!!
!!
!!
!!
!!
!!
!! 
!!
!!***


subroutine Grid_putPointData(blockid, cellType, variable, beginCount, &
     position, datablock)

  use Grid_data, ONLY : gr_iguard, gr_jguard, gr_kguard, gr_maxPatches


  implicit none
#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: blockid, variable, beginCount, cellType
  integer, dimension(MDIM), intent(in) :: position
  real, intent(in) :: datablock

  integer :: levelNum, patchNum

  integer, dimension(MDIM) :: posSamrai
  integer :: varIDSamrai

  !translate blockID to patchNum and levelNum
  levelNum = blockID / gr_maxPatches

  patchNum = mod(blockID, gr_maxPatches)


  !katie, use beginCount to make indices work for samrai
  !first interior cell is 0 in samrai
  if (beginCount == INTERIOR) then
     posSamrai(1) = position(1) - gr_iguard
     posSamrai(2) = position(2) - gr_jguard
     posSamrai(3) = position(3) - gr_kguard
  end if

  !do some error checking
  if (NDIM < 3) then
     posSamrai(3) = 1
  end if

  if (NDIM < 2) then
     posSamrai(2) = 1
  end if


  !samrai indexing starts at 0!
  varIDSamrai = variable - 1
  

  call samrai_put_point_data(patchNum, levelNum, varIDSamrai, posSamrai, datablock)




  
  return
end subroutine Grid_putPointData
