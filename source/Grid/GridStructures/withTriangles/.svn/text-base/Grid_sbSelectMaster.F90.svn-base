!!****if* source/Grid/GridStructures/withTriangles/Grid_sbSelectMaster
!!
!! NAME
!!  gr_sbSelectMaster
!!
!! SYNOPSIS
!!
!!  gr_sbSelectMaster()
!!  
!! DESCRIPTION 
!!  
!!  This routine is called from Grid_initDomain and from Driver_evoloveFlash
!!  for moving body. It identifies the processor that overlaps the most
!!  with the solid body as the master processor.
!!
!!
!! ARGUMENTS 
!!
!!***

#include "constants.h"
#include "Flash.h"

Subroutine Grid_sbSelectMaster()
  use Logfile_interface, ONLY : Logfile_open, Logfile_close
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkBoundBox
  use Grid_data, ONLY : gr_meshMe, gr_meshComm
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, gr_sbDebug, aelem
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Timers_interface, ONLY : Timers_start, Timers_stop
  implicit none
  include "Flash_mpi.h"

  real, dimension(2,MDIM) :: boundBox
  real, dimension(2) :: overlapTuple, maxOverlapTuple
  integer,dimension(MAXBLOCKS) :: listOfBlocks
  integer :: count, i, j, k, ierr
  integer :: logUnit, numCentroids, blkID
  logical, parameter :: localLogFile = .true.


  call Timers_start("body_select_master")
  call Grid_getListOfBlocks(LEAF, listOfBlocks, count)
  
  do k = 1, gr_sbNumBodies

     overlapTuple(1) = 0.0
     if (gr_sbBodyInfo(k) % comm /= MPI_COMM_NULL) then         

        numCentroids = size(gr_sbBodyInfo(k) % triangleCentroids,2)
        do j = 1, numCentroids
           gr_sbBodyInfo(k) % triangleCentroids(IAXIS,j) = &
             1./3. * (gr_sbBodyInfo(k) % xb(aelem(1,j)) + gr_sbBodyInfo(k) % xb(aelem(2,j)) + gr_sbBodyInfo(k) % xb(aelem(3,j)))

           gr_sbBodyInfo(k) % triangleCentroids(JAXIS,j) = &
                1./3. * (gr_sbBodyInfo(k) % yb(aelem(1,j)) + gr_sbBodyInfo(k) % yb(aelem(2,j)) + gr_sbBodyInfo(k) % yb(aelem(3,j)))
           
           gr_sbBodyInfo(k) % triangleCentroids(KAXIS,j) = &
                1./3. * (gr_sbBodyInfo(k) % zb(aelem(1,j)) + gr_sbBodyInfo(k) % zb(aelem(2,j)) + gr_sbBodyInfo(k) % zb(aelem(3,j)))
           do i = 1, count
              blkID = listOfBlocks(i)
              call Grid_getBlkBoundBox(blkID, boundBox)

              if (all(&
                   gr_sbBodyInfo(k) % triangleCentroids(1:NDIM,j) < boundBox(HIGH,1:NDIM) .and. &
                   gr_sbBodyInfo(k) % triangleCentroids(1:NDIM,j) > boundBox(LOW,1:NDIM))) then           
                 overlapTuple(1) = overlapTuple(1) + 1.0
              end if
           end do
        end do

        call MPI_Comm_rank(gr_sbBodyInfo(k) % comm, gr_sbBodyInfo(k) % myPE, ierr)
        overlapTuple(2) = gr_sbBodyInfo(k) % myPE
        call MPI_Comm_size(gr_sbBodyInfo(k) % comm, gr_sbBodyInfo(k) % numProcs, ierr)

        ! If only one process in the communicator
        if (gr_sbBodyInfo(k) % numProcs == 1) then
           gr_sbBodyInfo(k) % bodyMaster = gr_sbBodyInfo(k) % myPE
        else
           !The processor with blocks that overlap the most with the solid body
           !is the solid body master processor.
           call MPI_AllReduce(overlapTuple, maxOverlapTuple, 1, & 
                MPI_2Double_Precision, MPI_MAXLOC, gr_sbBodyInfo(k) % comm, ierr)
           gr_sbBodyInfo(k) % bodyMaster = int(maxOverlapTuple(2))
        endif
     end if

     if (gr_sbDebug) then
        call Logfile_open(logUnit,localLogFile)
        write(logUnit,'(a,i8,a,es14.6,a,l8)') "Body", k, &
             ", overlap", overlapTuple(1), &
             ", in communicator", &
             gr_sbBodyInfo(k) % comm /= MPI_COMM_NULL
        call Logfile_close(localLogFile)
     end if

  end do

  call Timers_stop("body_select_master")

End Subroutine Grid_sbSelectMaster
