!!****f* source/Grid/Grid_getBlkIDFromPos
!!
!! NAME
!!  Grid_getBlkIDFromPos
!!
!! SYNOPSIS
!!
!!  call Grid_getBlkIDFromPos(real(IN)     :: pos(:), 
!!                            integer(IN) :: blkList(:),
!!                            integer(IN) :: blkCount,
!!                            integer(OUT) :: blkID,
!!                            integer(OUT) :: procID,
!!                   optional,integer(IN)  :: comm)
!!  
!! DESCRIPTION 
!! 
!!  Returns the processor and blkID
!!  containing the cell that overlaps with the 
!!  specified position co-ordinate
!! 
!! 
!! ARGUMENTS 
!!
!!  pos        :: co-ordinates of the point
!!  blkList    :: the list of blocks to search
!!  blkCount   :: the count of blocks in the list
!!  blkID      :: the local blockID of the block that contains the point
!!  procID     :: the processor ID that contains the point
!!  comm       :: if a communicator other than the default mesh communicator is
!!                desired, it should be specified here
!!
!! EXAMPLE
!!
!!***


#include "constants.h"

subroutine Grid_getBlkIDFromPos(pos,blkList, blkCount,blockId, procID,comm)

  implicit none

  real, dimension(1:MDIM), intent(IN) :: pos
  integer,intent(IN)  :: blkCount
  integer,dimension(blkCount),intent(IN) :: blkList
  integer, intent(OUT) :: blockID, procID
  integer,optional,intent(IN) :: comm
  
  procID=NONEXISTENT
  blkID=NONEXISTENT
  
end subroutine Grid_getBlkIDFromPos
