!!****if* source/Grid/GridStructures/withTriangles/pointTopoint/bitmap/unitTest/Grid_solidBodyUnitTest
!!
!! NAME
!!
!!  Grid_solidBodyUnitTest
!!
!! SYNOPSIS
!!
!!  Grid_solidBodyUnitTest(integer, intent(in):: fileUnit,
!!                logical, intent(inout)::perfect  )
!!
!! DESCRIPTION
!!
!!  This routine tests the particle exchange information. Specifically it
!!  checks that the proc ID property of each particle matches the processor 
!!  it belongs to. The proc ID property of each particle is sent from its 
!!  respective processor to the master processor.
!!
!! ARGUMENTS
!!  fileUnit - open f90 write unit
!!  perfect - returns a true if the test passed, false otherwise
!! NOTES
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine Grid_solidBodyUnitTest(fileUnit,perfect)
  use Grid_interface, ONLY : Grid_getBlkBoundBox, Grid_getListOfBlocks
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, totalPart
  use Grid_data, ONLY : gr_meshComm
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  include "Flash_mpi.h"
  integer, intent(IN)           :: fileUnit ! Output to file
  logical, intent(INOUT)        :: perfect  ! Flag to indicate errors

  real, dimension(2,MDIM) :: boundBox
  integer, dimension(MAXBLOCKS) :: listOfBlocks
  real, dimension(MDIM) :: pos
  integer :: p, count, b, blkID, proc, i, ierr, myPE, numProcs
  logical :: expectedParticleProcID
  logical :: perfectLocal

  call Timers_start("body_check_results")
  perfectLocal = .true.

  call Grid_getListOfBlocks(LEAF, listOfBlocks, count)

  !Loop over each created particle.  If the particle is within
  !a block on the master processor then the proc ID property
  !should be equal to the master processor, otherwise it should
  !not be equal to the master processor and also be positive.
  do b = 1, gr_sbNumBodies
     if (gr_sbBodyInfo(b) % myPE == gr_sbBodyInfo(b) % bodyMaster) then
        totalPart = gr_sbBodyInfo(b) % totalPart

        call MPI_Comm_rank(gr_meshComm, myPE, ierr)
        call MPI_Comm_size(gr_meshComm, numProcs, ierr)
        if (myPE /= gr_sbBodyInfo(b) % myPE) then
           call Driver_abortFlash("Process ID mistake")
        end if

        do p = 1, totalPart
           pos(IAXIS) = gr_sbBodyInfo(b) % particles(POSX_PART_PROP,p)
           pos(JAXIS) = gr_sbBodyInfo(b) % particles(POSY_PART_PROP,p)
           pos(KAXIS) = gr_sbBodyInfo(b) % particles(POSZ_PART_PROP,p)
           proc = gr_sbBodyInfo(b) % particles(PROC_PART_PROP,p)
           
           expectedParticleProcID = .false.
           do i = 1, count
              blkID = listOfBlocks(i)
              call Grid_getBlkBoundBox(blkID, boundBox)

!                 print *, "PROC:", gr_sbBodyInfo(b) % myPE, "block:", boundbox(1:2,1:2)
                 !If within a block on the master processor then we
                 !expect procID to be the master processor.
              if (all(&
                   pos(1:NDIM) >= boundBox(LOW,1:NDIM) .and. &
                   pos(1:NDIM) < boundBox(HIGH,1:NDIM))) then
                 if (proc == gr_sbBodyInfo(b) % myPE) then
                    expectedParticleProcID = .true.
                    exit
                 end if
              end if
           end do

              !If not on a block on the master processor then we
              !expect a different processor ID.
           if (.not.expectedParticleProcID) then
              if (proc /= gr_sbBodyInfo(b) % myPE .and. &
                   proc >= 0 .and. proc < numProcs) then
                 expectedParticleProcID = .true.
              end if
           end if
           
           if (.not.expectedParticleProcID) then
!                 write(fileUnit,*) 'particle p = ', p, &
!                      ' has unexpected processor ID of ', proc
              write(fileUnit,*) 'Body no = ', b, 'particle p = ', p, &
                   ' has unexpected processor ID of ', proc, 'pos =', pos
              perfectLocal = .false.
           end if
        end do
     end if
  end do

  call MPI_AllReduce(perfectLocal, perfect, 1, MPI_LOGICAL, MPI_LAND, &
       gr_meshComm, ierr)

  call Timers_stop("body_check_results")

end subroutine Grid_solidBodyUnitTest
