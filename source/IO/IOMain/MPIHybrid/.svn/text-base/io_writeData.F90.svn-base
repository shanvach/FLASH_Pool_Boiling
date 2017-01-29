!!****if* source/IO/IOMain/MPIHybrid/io_writeData
!!
!! NAME
!!
!!  io_writeData
!!
!! SYNOPSIS
!!
!!  io_writeData(integer,intent(IN)  :: fileID)
!!
!! DESCRIPTION
!!
!!  This function writes the checkpoint data to a file to store the 
!!  state of the simulation. It uses MPI_IO, and outputs the data
!!  in binary format. Users should exercise caution when using this
!!  data since it is "endian" sensitive.
!!
!! ARGUMENTS
!!
!!   fileID : file handle
!!
!!
!!
!!***


subroutine io_writeData(fileID)
  use IO_data, ONLY : io_dataType
  use Grid_interface, ONLY : Grid_getBlkIndexLimits
  use physicalData, ONLY : unk
implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  
  integer,intent(IN) :: fileID
  integer :: ioTypeSize, flashReal_typeSize
  integer :: count
  integer(kind=8)::offset
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: blockID=1

  call MPI_Type_size(io_dataType,ioTypeSize,ierr)
  call MPI_Type_size(FLASH_REAL,flashReal_typeSize,ierr)
  
  count = ioTypeSize/flashReal_typeSize

  offset = 0
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

  call MPI_File_set_view(fileID,MPI_BYTE,io_dataType,"native",&
       MPI_INFO_NULL,ierr)

  do i = UNK_VARS_BEGIN,UNK_VARS_END
     call MPI_File_write_all&
          (fileID,unk(i,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                        blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
                        blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)),&
                        count,FLASH_REAL,status,ierr)
  end do

  return
end subroutine Grid_dump
