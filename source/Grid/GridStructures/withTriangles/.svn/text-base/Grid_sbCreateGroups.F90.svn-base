!!****if* source/Grid/GridStructures/withTriangles/Grid_sbCreateGroups
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
!!  This routine is called from Grid_initDomain and from Drivier_evolveFlash
!!  for moving body. It creates communicators
!!  for each solid body that overlap with it.
!!
!!
!! ARGUMENTS 
!!
!!***

#include "constants.h"
#include "Flash.h"

Subroutine Grid_sbCreateGroups()
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
     gr_sbBodyInfo(j) % myPE = NONEXISTENT
     gr_sbBodyInfo(j) % numProcs = NONEXISTENT
     gr_sbBodyInfo(j) % bodyMaster = NONEXISTENT
  end do

  call Timers_stop("body_create_comms")

End Subroutine Grid_sbCreateGroups
