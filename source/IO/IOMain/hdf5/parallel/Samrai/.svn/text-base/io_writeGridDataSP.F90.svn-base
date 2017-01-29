!!****if* source/IO/IOMain/hdf5/parallel/Samrai/io_writeGridDataSP
!!
!! NAME
!!
!!  io_writeGridDataSP
!!
!! 
!! SYNOPSIS
!!
!!  io_writeGridDataSP(filenum, simtime)
!!
!!  io_writeGridDataSP(integer, real)
!!
!!
!! DESCRIPTION
!!
!!  Dump out a plotfile using parallel HDF5.  The IO is done in parallel -- 
!!  no copying of the data to a single processor to do the writing is 
!!  performed.  
!!
!!  This is the SINGLE PRECISION version of the plotfile -- temporary
!!  storage is used to recast a variable (for every zone/block) into
!!  single precision before passing it onto the SP version of the C HDF 5
!!  write routines.
!! 
!!  The data for all blocks is recast and written together.  This makes the
!!  amount of data that is written very large, which should perform better
!!  on the parallel filesystems.  The overhead for storing an entire 
!!  variable (with corners) is small, <~ 1%.
!!
!!  HDF5 uses MPI-IO (via ROMIO) to support parallel IO.  Each processor
!!  must open the file, create the datasets, and dataspaces for each HDF
!!  record.  
!!
!!  A single record for each of the PARAMESH data structures is created.  A
!!  processor only writes to a subset of this record.  Each record has a 
!!  dimension with length = tot_blocks.  The offset of a processor into this 
!!  dimension is computed by looking at the total number of blocks that are
!!  below the current processor.
!!
!!  In this version of the plotfile, each variable is given its own 
!!  record -- this makes it easier to change the variable list in the
!!  future without disturbing the format of the file.
!!
!!  The include file -- hdf5_flash.h is used for the C routines and mirrors
!!  the necessary data from physicaldata.fh
!!
!!
!! NOTES
!!
!!  This version of the plotfile routine requires HDF5 v. 1.4.0 or later.
!!
!!
!! ARGUMENTS
!!
!!  filenum        The number of the file to output
!!
!!  simtime        The current simulation time
!!
!!***

subroutine io_writeGridDataSP (myPE, numProcs, fileID)

implicit none
  integer, intent(in) :: myPE, numProcs, fileID

  return
end subroutine io_writeGridDataSP

