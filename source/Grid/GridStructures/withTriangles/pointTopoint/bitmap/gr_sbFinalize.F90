!!****if* source/Grid/GridStructures/withTriangles/pointTopoint/bitmap/gr_sbFinalize
!!
!! NAME
!!
!!  gr_sbFinalize
!!
!! SYNOPSIS
!!
!!  gr_sbFinalize()
!!
!! DESCRIPTION
!!
!!  Called from Grid_finalize. Deallocates the data structure that holds all information
!!  about the solid bodies.
!!
!! ARGUMENTS
!!
!!
!!***


Subroutine gr_sbFinalize()
use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, gr_sbParticleCount
implicit none
include "Flash_mpi.h"
integer :: b, ierr

do b = 1, gr_sbNumBodies
   if (gr_sbBodyInfo(b) % myPE == gr_sbBodyInfo(b) % bodyMaster) then
         deallocate(gr_sbBodyInfo(b) % particles)      
   else
      if(gr_sbParticleCount(b) > 0) then
         deallocate(gr_sbBodyInfo(b) % particles)
      end if
   endif
end do
if (allocated(gr_sbParticleCount)) deallocate(gr_sbParticleCount)
!deallocate(gr_sbBodyInfo)

End Subroutine gr_sbFinalize
