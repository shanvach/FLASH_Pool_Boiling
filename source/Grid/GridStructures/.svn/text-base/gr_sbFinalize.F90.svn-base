!!****if* source/Grid/GridStructures/gr_sbFinalize
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
use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies
implicit none
include "Flash_mpi.h"
integer :: b, ierr

do b = 1, gr_sbNumBodies
   if (gr_sbBodyInfo(b) % comm /= MPI_COMM_NULL) then
      call MPI_Comm_free(gr_sbBodyInfo(b) % comm, ierr)
      deallocate(gr_sbBodyInfo(b) % particles)
   end if
end do

deallocate(gr_sbBodyInfo)

End Subroutine gr_sbFinalize
