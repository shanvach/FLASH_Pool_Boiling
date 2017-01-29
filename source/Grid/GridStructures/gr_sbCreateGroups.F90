!!****if* source/Grid/GridStructures/gr_sbCreateGroups
!!
!! NAME
!!  gr_sbCreateGroups
!!
!! SYNOPSIS
!!
!!  gr_sbCreateGroups()
!!  
!! DESCRIPTION 
!!  
!!  This routine is called from Grid_initDomain. It creates communicators
!!  for each solid body and identifies the processor that overlaps the most
!!  with the solid body as the master processor.
!!
!!
!! ARGUMENTS 
!!
!!***

#include "constants.h"
#include "Flash.h"

Subroutine gr_sbCreateGroups()
  use Logfile_interface, ONLY : Logfile_open, Logfile_close
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkBoundBox
  use Grid_data, ONLY : gr_meshMe, gr_meshComm
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, gr_sbDebug
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Timers_interface, ONLY : Timers_start, Timers_stop
  implicit none
  include "Flash_mpi.h"

  real, dimension(2,MDIM) :: boundBox
  real, dimension(MDIM) :: lowOverlap, highOverlap
  real, dimension(2) :: overlapTuple, maxOverlapTuple
  integer,dimension(MAXBLOCKS) :: listOfBlocks
  integer :: count, i, j, blkID, ierr
  integer :: color, key, logUnit
  logical, parameter :: localLogFile = .true.


  call Timers_start("body_create_comms")
  call Grid_getListOfBlocks(LEAF, listOfBlocks, count)
  
  do j = 1, gr_sbNumBodies
     overlapTuple(1) = 0.0

     do i = 1, count
        blkID = listOfBlocks(i)
        call Grid_getBlkBoundBox(blkID, boundBox)

        if (all(&
             gr_sbBodyInfo(j) % boundBox(LOW,1:NDIM) < boundBox(HIGH,1:NDIM) .and. &
             gr_sbBodyInfo(j) % boundBox(HIGH,1:NDIM) > boundBox(LOW,1:NDIM))) then
           
           !This will not work for periodic boundary conditions.
           lowOverlap(1:NDIM) = max(&
                gr_sbBodyInfo(j) % boundBox(LOW,1:NDIM),&
                boundBox(LOW,1:NDIM))

           highOverlap(1:NDIM) = min(&
                gr_sbBodyInfo(j) % boundBox(HIGH,1:NDIM),&
                boundBox(HIGH,1:NDIM))
           
           overlapTuple(1) = overlapTuple(1) + &
                product(highOverlap(1:NDIM) - lowOverlap(1:NDIM))

        end if
     end do

     key = gr_meshMe
     if (overlapTuple(1) == 0.0) then
        !The MPI_Comm_split call will return MPI_COMM_NULL communicator
        color = MPI_UNDEFINED
     else
        color = 1
     end if
     call MPI_Comm_split(gr_meshComm, color, key, gr_sbBodyInfo(j) % comm, ierr)


     nullify(gr_sbBodyInfo(j) % particles)
     if (gr_sbBodyInfo(j) % comm == MPI_COMM_NULL) then 
        gr_sbBodyInfo(j) % myPE = NONEXISTENT
        gr_sbBodyInfo(j) % numProcs = NONEXISTENT
        gr_sbBodyInfo(j) % bodyMaster = NONEXISTENT
     else
        call MPI_Comm_rank(gr_sbBodyInfo(j) % comm, gr_sbBodyInfo(j) % myPE, ierr)
        overlapTuple(2) = gr_sbBodyInfo(j) % myPE
        call MPI_Comm_size(gr_sbBodyInfo(j) % comm, gr_sbBodyInfo(j) % numProcs, ierr)

        ! If only one process in the communicator
        if (gr_sbBodyInfo(j) % numProcs == 1) then
           gr_sbBodyInfo(j) % bodyMaster = gr_sbBodyInfo(j) % myPE
        else        
           !The processor with blocks that overlap the most with the solid body
           !is the solid body master processor.
           call MPI_AllReduce(overlapTuple, maxOverlapTuple, 1, & 
                MPI_2Double_Precision, MPI_MAXLOC, gr_sbBodyInfo(j) % comm, ierr)
           gr_sbBodyInfo(j) % bodyMaster = int(maxOverlapTuple(2))

           if (gr_sbBodyInfo(j) % myPE == gr_sbBodyInfo(j) % bodyMaster) then
           end if
        endif
     end if

     if (gr_sbDebug) then
        call Logfile_open(logUnit,localLogFile)
        write(logUnit,'(a,i8,a,es14.6,a,l8)') "Body", j, &
             ", overlap", overlapTuple(1), &
             ", in communicator", &
             gr_sbBodyInfo(j) % comm /= MPI_COMM_NULL
        call Logfile_close(localLogFile)
     end if

  end do

  call Timers_stop("body_create_comms")

End Subroutine gr_sbCreateGroups
