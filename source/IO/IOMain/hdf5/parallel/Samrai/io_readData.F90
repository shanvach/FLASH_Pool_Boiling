!!****if* source/IO/IOMain/hdf5/parallel/Samrai/io_readData
!!
!! NAME
!!
!!  io_readData
!!
!!
!! SYNOPSIS
!!
!!  io_readData(io_globalMe, numProcs)
!!
!!  io_readData(integer, integer)
!!
!!
!! DESCRIPTION
!!
!!  This is the reading counterpart to io_writeData.  It reads an HDF5
!!  file and distributes it to the processors to restart a simulation.
!!
!!  All reading is done using parallel HDF5 -- no explicit data movement
!!  is used.
!!
!!***




subroutine io_readData(io_globalMe, numProcs)

implicit none
  integer, intent(in) :: io_globalMe, numProcs

  return
end subroutine io_readData
