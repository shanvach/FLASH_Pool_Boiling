!!****if* source/IO/IOMain/MeshDump/IO_writeCheckpoint
!!
!! NAME
!!
!!  IO_writeCheckpoint
!!
!!
!! SYNOPSIS
!!
!!  IO_writeCheckpoint()
!!
!!
!!
!! DESCRIPTION
!!
!!  This is a generic call to write the important simulation data to a
!!  checkpoint file.  A checkpoint file writes a few different types of
!!  data to a file, first the physical data like, temperature, pressure, density
!!  etc. which are stored in each cell on the grid.  
!!  Second, in order to recreate the simulation from a checkpoint file a
!!  number of other single quantities are needed as well.  We call these
!!  scalar values which include simTime, dt, nstep, globalNumBlocks etc.
!!  We also store descriptive strings that describe the simulation run.
!!
!!  The same IO_writeCheckpoint routine is called regardless of the type of
!!  file being written, (such as hdf5 parallel, hdf5 serial or pnetcdf)
!!  IO_writeCheckpoint prepares the Grid_ioData (like getting the
!!  globalNumBlocks) and collects the scalars wanting to be checkpointed
!!  from each unit. 
!!  IO_writeCheckpoint then calls four methods, io_initFile, io_writeData, 
!!  IO_writeParticles and io_closeFile.  Each of these routines _is_ specific
!!  to the type of io library used and have their own implementation.  
!!  In addition, io_writeData has its own
!!  implementation for io library and type of grid (UG, Paramesh, or other)
!!
!!
!!  In FLASH IO_writeCheckpoint is called from IO_output (or IO_outputInitial or
!!  IO_outputFinal) IO_output checks to see if enough wall clock time,
!!  simTim, or nsteps has passed to checkpoint.
!!
!!  We have put IO_writeCheckpoint in the API because a user may want to write
!!  a checkpoint at another time or for another reason without having to go through
!!  IO_output.  For most flash users IO_writeCheckpoint will only ever be
!!  called through IO_output
!!
!! ARGUMENTS
!! 
!!
!! NOTES
!!
!!  For those familiar with FLASH2, breaking up the checkpoint routine into
!!  these four different methods is a change.  Because FLASH3 now supports
!!  different grid packages and we are committed to supporting both
!!  hdf5 and parallel netCDF having each grid and io library writing its
!!  own checkpoint file proved to be a lot of code duplication.  We believe
!!  that while dividing up the checkpoint routines created more files it 
!!  will in the end be easier to maintain.
!!
!!***
subroutine IO_writeCheckpoint()
  use Logfile_interface, ONLY : Logfile_stamp
  use IO_data, ONLY : io_globalMe, io_globalNumProcs
  implicit none

#include "constants.h"  
#include "Flash.h"  

  ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(1,2) :: strBuff
  character(len=MAX_STRING_LENGTH) :: filenamePrefix
  integer, save :: chkNum = 0

  interface
     subroutine io_meshDump(myPE, numProcs, filenamePrefix)
       implicit none
       integer, intent(IN) :: myPE, numProcs
       character(len=MAX_STRING_LENGTH), intent(IN) :: filenamePrefix
     end subroutine io_meshDump
  end interface
  
  write(filenamePrefix, fmt='(''CENTER'',i3.3,''_'')') chkNum
  call io_meshDump(io_globalMe, io_globalNumProcs, filenamePrefix)  !Add CENTER arg plus variable list

  if (io_globalMe == MASTER_PE) then
     write (strBuff(1,1), "(A)") "type"
     write (strBuff(1,2), "(A)") "checkpoint"
     call Logfile_stamp(strBuff, 1, 2, "[IO_writeCheckpoint] close")
     print *, '*** Wrote checkpoint file ***' 
  end if

  chkNum = chkNum + 1
end subroutine IO_writeCheckpoint
