!!****if* source/Simulation/SimulationMain/Flash2Convert/IO_outputInitial
!!
!! NAME
!!
!!  IO_outputInitial
!!
!!
!! SYNOPSIS
!!
!!
!!  IO_outputInitial() 
!!                   
!!                   integer(in) :: nbegin,
!!                   real(in) :: initialSimTime  
!!                  
!!
!!
!! DESCRIPTION
!!
!!  Outputs the converted Flash 3 file as the same checkpoint as is being read 
!!  in from Flash 2 of a Flash 3 run.
!!  Also removed the call to output a plot file (does not make sense in this 
!!  case).
!!
!!
!!  This routine is called before the main timestep loop.  It outputs the 
!!  initial data to a checkpoint file and plotfile, and particle plotfiles
!!
!!  If particles are not included a stub (empty) routine will be called.
!!
!!
!! ARGUMENTS
!!
!!  myPE - current processor
!!  numProcs - number of processors running the simulation
!!  nbegin - initial step of simulation
!!  initialSimTime - initial simulation time
!!
!!
!!***


subroutine IO_outputInitial(myPE, numProcs, nbegin, initialSimTime)

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use IO_data, ONLY : io_integralFreq, io_lastWallClockCheckpoint, &
       io_redshift, io_justCheckpointed, io_restart, io_checkpointFileNumber
  use IO_interface, ONLY : IO_writeIntegralQuantities, &
    IO_writeCheckpoint, IO_writePlotfile, IO_writeParticles

  implicit none

#include "Flash_mpi.h"

  integer, intent(in) :: myPE, numProcs, nbegin
  real, intent(in) :: initialSimTime

  !write the diagnostic quantities for the .dat file
  if(io_integralFreq > 0) then
     call IO_writeIntegralQuantities(myPE, 1, initialSimTime)
  end if



  if(io_restart) then
     call RuntimeParameters_get('checkpointFileNumber', io_checkpointFileNumber)
     call IO_writeCheckpoint()
     io_lastWallClockCheckpoint = MPI_Wtime()
     io_justCheckpointed = .true.
  end if

  !call IO_writePlotfile()

  !call IO_writeParticles( .false.)

end subroutine IO_outputInitial
