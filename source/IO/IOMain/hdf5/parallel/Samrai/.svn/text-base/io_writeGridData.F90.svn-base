!!****if* source/IO/IOMain/hdf5/parallel/Samrai/io_writeGridData
!!
!! NAME
!!
!!  io_writeGridData
!!
!!
!! SYNOPSIS
!!
!!  io_writeGridData(int myPE, int numProcs, int fileID)
!!           
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
!!  dimension with length = tot_blocks.  The offset of a processor into this
!!  dimension is computed by looking at the total number of blocks that are
!!  below the current processor.
!!
!!  In this version of io_writeData, each variable is given its own
!!  record -- this makes it easier to change the variable list in the
!!  future without disturbing the format of the file.
!!
!!  
!!
!!  
!!
!! ARGUMENTS
!! 
!!  myPE - current processor number
!!  numProcs - number of procs running current simulation
!!  fileID - integer file identifier for hdf5 file
!!  nzones_block - holds nxb, nyb, nzb
!!
!!***


subroutine io_writeGridData(myPE, numProcs, fileID) 
  
implicit none
  integer, intent(in) :: myPE, numProcs, fileID

end subroutine io_writeGridData
