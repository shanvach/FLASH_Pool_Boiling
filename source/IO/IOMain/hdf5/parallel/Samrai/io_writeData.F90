!!****if* source/IO/IOMain/hdf5/parallel/Samrai/io_writeData
!!
!! NAME
!!
!!  io_writeData
!!
!!
!! SYNOPSIS
!!
!!  io_writeData() 
!!              
!!             integer(in) :: fileID) 
!!          
!!
!!
!!
!! DESCRIPTION
!!
!!  This function writes the checkpoint data to an hdf5 file to store the 
!!  paramesh data.  IO is done in parallel -- no copying of the data to 
!!  a single processor
!!  to do the writing is performed.  HDF5 v. 1.4.0 or later is required
!!
!!  HDF5 uses MPI-IO (via ROMIO) to support parallel IO.  Each processor
!!  must open the file, create the datasets, and dataspaces for each HDF
!!  record.
!!
!!  A single record for each of the PARAMESH data structures is created.  A
!!  processor only writes to a subset of this record.  Each record has a
!!  dimension with length = globalNumBlocks.  The offset of a processor into this
!!  dimension is computed by looking at the total number of blocks that are
!!  below the current processor.
!!
!!  In this version of io_writeData, each variable is given its own
!!  record -- this makes it easier to change the variable list in the
!!  future without disturbing the format of the file.
!!
!!
!! ARGUMENTS
!! 
!!  io_globalMe - current processor number
!!  numProcs - number of procs running current simulation
!!  fileID - integer file identifier for hdf5 file
!!
!!
!! NOTES
!!  variables that start with "io_" belong to the data module IO_data.
!!  the "io_" is meant to indicate that these variables have IO unit 
!!  scope.   For performance purposes IO particularly, io_writeData
!!  and read data are allowed to access unk directly.  Additionally,
!!  io_writeData needs access to some variables that are specific to 
!!  Grid.  These are stored in Grid_ioData and variables associated
!!  with Grid_ioData start with "gio_"
!!
!!***


subroutine io_writeData(io_globalMe, numProcs, fileID) 

implicit none
  integer, intent(in) :: myPe, numProcs, fileID
  

  return

end subroutine io_writeData
