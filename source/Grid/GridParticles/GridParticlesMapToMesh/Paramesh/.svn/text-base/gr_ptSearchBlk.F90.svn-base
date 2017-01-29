!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/Paramesh/gr_ptSearchBlk
!!
!! NAME
!!  gr_ptSearchBlk
!!
!! SYNOPSIS
!!
!!  gr_ptSearchBlk(integer,dimension(MDIM), intent(IN) :: cornerID, &
!!                 integer,dimension(BLKID:REFLEVELDIF), intent(INOUT) :: negh)
!!
!! DESCRIPTION
!!
!!  Routine to check whether a paricular block on this processor matches the 
!!  corner ID argument.  If there is a match, the negh array is updated.
!!
!! ARGUMENTS
!!               cornerID:   The corner ID we will try to match.
!!               negh:   Array containing information about the destination block(s)
!!
!! PARAMETERS
!! 
!!***


subroutine gr_ptSearchBlk(cornerID,negh)

#include "Flash.h"
#include "constants.h"
#include "gr_ptMapToMesh.h"

#ifdef FLASH_GRID_PARAMESH
  use tree, ONLY : lnblocks,nodetype

  use Grid_data, ONLY : gr_oneBlock, gr_meshMe
#endif
  implicit none  

  integer,dimension(MDIM), intent(IN) :: cornerID
  integer,dimension(BLKID:REFLEVELDIF), intent(INOUT) :: negh
  integer :: lb
  logical :: notFound
  integer, dimension(NDIM) :: diff
  
#ifdef FLASH_GRID_PARAMESH  
  lb=0
  
  notFound=.true.
  do while ((lb<lnblocks).and.notFound)
     lb=lb+1
     if(nodetype(lb)==LEAF) then
        diff=cornerID(1:NDIM)-gr_oneBlock(lb)%cornerID(1:NDIM)
        if(maxval(abs(diff))==0)then
           negh(BLKID)=lb
           negh(BLKPROC)=gr_meshMe
           notFound=.false.
        end if
     end if
  end do

  return
#endif
end subroutine gr_ptSearchBlk
