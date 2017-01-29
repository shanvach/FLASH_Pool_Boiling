!!****if* source/Simulation/SimulationMain/Flash2Convert/Driver_evolveFlash
!!
!! NAME
!!
!!  Driver_evolveFlash
!!
!! SYNOPSIS
!!
!!  Driver_evolveFlash()
!!
!! DESCRIPTION
!!
!! This version of Driver_evolveFlash is meant for use only with the Flash 2 to 
!! Flash 3 file converter.  It may be usable for other setups that require that 
!! Flash 3 not process beyond initialization.  This is mostly a do-nothing stub,
!! that notes in the logfile that the function has been called (this may help 
!! in performing diagnostics).
!!
!!  
!!
!! NOTES
!!
!! The Driver unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. These variables begin with "dr_"
!! like, dr_globalMe or dr_dt, dr_beginStep, and are stored in fortran
!! module Driver_data (in file Driver_data.F90. The other variables
!! are local to the specific routine and do not have the prefix "dr_"
!!
!!
!!***


#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif


subroutine Driver_evolveFlash()

  use Driver_data, ONLY: dr_globalMe, dr_numProcs, dr_nbegin, &
       dr_nend, dr_dt, dr_wallClockTimeLimit, &
       dr_tmax, dr_simTime, dr_redshift, &
       dr_nstep, dr_dtOld, dr_dtNew, dr_restart, dr_elapsedWCTime
  use Driver_interface, ONLY : Driver_sourceTerms, Driver_computeDt
  use Logfile_interface, ONLY : Logfile_stamp, Logfile_close
  use Timers_interface, ONLY : Timers_start, Timers_stop, &
    Timers_getSummary
  use Particles_interface, ONLY : Particles_advance, Particles_dump
  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
    Grid_getListOfBlocks, Grid_updateRefinement
  use Hydro_interface, ONLY : Hydro
  use Gravity_interface, ONLY :  Gravity_potentialListOfBlocks
  use IO_interface, ONLY :IO_output,IO_outputFinal

  implicit none

#include "constants.h"
#include "Flash.h"
 include "Flash_mpi.h"

  integer   :: localNumBlocks

  integer :: blockCount, ierr
  integer :: blockList(MAXBLOCKS)

  
  ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(3,2) :: strBuff
  character(len=15) :: numToStr

  call Logfile_stamp( 'Entering evolution loop' , '[Driver_evolveFlash]')


  call Timers_start("evolution")


  call Timers_stop("evolution")

  call Logfile_stamp( 'Exiting evolution loop' , '[Driver_evolveFlash]')


  call Timers_getSummary(dr_globalMe, dr_nstep)


  call Logfile_stamp( "FLASH run complete.", "LOGFILE_END")

  call Logfile_close()


  return
  
end subroutine Driver_evolveFlash



