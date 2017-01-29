!!****if* source/IO/IOMain/MPIHybrid/io_getOutPutName
!!
!! NAME
!!  io_getOutPutName
!!
!! SYNOPSIS
!!
!!  io_getOutPutName(integer, intent(IN) ::filenum, 
!!                   character(:), intent(IN) :: strname, 
!!                   character(:), intent(OUT) :: filename)
!!
!! DESCRIPTION
!!
!!  gets the name of the output file for a checkpoint or plotfile
!!  it takes an integer filenum along with a string describing the
!!  file, ie hdf5_chk_ or ncmpi_chk_ and combines it with the number in
!!  a string format 
!!
!! ARGUMENTS
!!
!!  filenum - integer representing the checkpoint or plotfile number
!!  strname - string descriptor of file like, hdf5_chk or hdf5_plt
!!  filename - string name given to the file
!!
!!***


subroutine io_getOutPutName(filenum, filetypestr, strname, filename)

  use IO_data, ONLY : io_baseName, io_outputDir, io_group

  implicit none
  
#include "constants.h"


  integer, intent(in) :: filenum
  character (len=MAX_STRING_LENGTH), intent(inout) :: filename
  character(len=*), intent(in) :: strname, filetypestr

  ! create a character variable to hold the string representation of the block
  ! number.  Note this is set to be 4 characters long (i.e. max = 9999).  
  character (len=4) ::  fnumString
  character (len=4) :: fgrpString
  integer :: pos



  ! if io_outputDir isn't empty, and it doesn't have a 
  ! directory seperator at the end, add the directory seperator.
  pos = index(io_outputDir, ' ')
  if (pos.gt.1) then
     if (io_outputDir(pos-1:pos-1).ne."/") then
        io_outputDir(pos:pos) = "/"
     end if
  end if
  
  ! if io_outputDir is current directory specified 
  ! with a './', just get rid of it.
  if (io_outputDir == "./") then
     io_outputDir = ""
  end if

  filename = " "
  write (fnumString, '(i4.4)') filenum
  write(fgrpString, '(i4.4)')io_group
  pos = index(io_baseName, ' ')
  filename = io_baseName(:pos-1) // filetypestr // strname  // fnumString &
                                // fgrpString
  pos = index(io_outputDir, ' ')
  filename = io_outputDir(:pos-1) // filename 


end subroutine io_getOutPutName
