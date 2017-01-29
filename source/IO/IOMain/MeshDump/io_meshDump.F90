!!****if* source/IO/IOMain/MeshDump/io_meshDump
!!
!! NAME
!!
!!  io_meshDump
!!
!!
!! SYNOPSIS
!!
!!  io_meshDump
!!
!! DESCRIPTION
!!
!!  This subroutine dumps mesh data to file for diagnostic purposes.
!!  It is designed so that we can verify FLASH UG gives identical
!!  results as Chombo UG.  A separate file is written for each processor.
!!
!! ARGUMENTS
!! 
!!  myPE - current processor number
!!  numProcs - number of procs running current simulation
!!
!!***

!!REORDER(4): dataPtr
#include "constants.h"
#include "Flash.h"

subroutine io_meshDump(myPE, nproc, fileNamePrefix)
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkIndexLimits, &
       Grid_getBlkPtr, Grid_releaseBlkPtr
  implicit none
  include 'Flash_mpi.h'
  integer, intent(IN) :: myPE, nproc
  character(len=MAX_STRING_LENGTH), intent(IN) :: fileNamePrefix

  real, dimension(:,:,:,:), pointer :: dataPtr
  integer,dimension(MAXBLOCKS) :: listOfBlocks
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: v, i, j, k, ierr, ii, blockID, globalBlockID, count
  integer, allocatable, dimension(:) :: blocks
  integer, parameter :: unit_txt = 12, unit_bin = 13
  character(len=MAX_STRING_LENGTH) :: fileNameTxt, fileNameBin, blkProcSuffix


  call Grid_getListOfBlocks(ALL_BLKS,listOfBlocks,count)

  allocate(blocks(nproc))
  blocks = 0
  blocks(myPE+1) = count
  call MPI_Allreduce(MPI_IN_PLACE, blocks, nproc, &
       FLASH_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

  do ii = 1, count
     blockID = listOfBlocks(ii)

     !Create a zero-based global block ID.
     if (myPE == 0) then
        globalBlockID = blockID - 1
     else
        globalBlockID = sum(blocks(1:myPE)) + blockID - 1
     end if

     !Open a unique file for each block by using the global block ID.
     write(blkProcSuffix, fmt='(''B'',i3.3,''_P'',i3.3)') globalBlockID, myPE
     fileNameTxt = trim(fileNamePrefix) // trim(blkProcSuffix) // ".txt"
     fileNameBin = trim(fileNamePrefix) // trim(blkProcSuffix) // ".bin"
     open(unit=unit_txt, file=fileNameTxt)
     open(unit=unit_bin, file=fileNameBin, form='unformatted')

     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
     call Grid_getBlkPtr(blockID, dataPtr)
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              do v = 1, NUNK_VARS
                 write(unit_txt,'(4i4, D20.12)') v, i, j, k, dataPtr(v,i,j,k)
                 write(unit_bin) dataPtr(v,i,j,k)
              end do
           end do
        end do
     end do
     call Grid_releaseBlkPtr(blockID, dataPtr)

     close(unit_txt)
     close(unit_bin)
  end do

  deallocate(blocks)
end subroutine io_meshDump
