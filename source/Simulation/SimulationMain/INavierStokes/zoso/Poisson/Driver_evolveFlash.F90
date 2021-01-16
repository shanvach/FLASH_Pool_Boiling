!!****if* source/Simulation/SimulationMain/INavierStokes/2D/bhagaWeber_mcHYPRE/Driver_evolveFlash
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
!!  This is the main global driver for simulations that are:
!!      Spatially refined, State form, strang split
!!
!!  DOC: Driver_evolveFlash needs more explanation 
!!
!! NOTES
!!
!!  variables that begin with "dr_" like, dr_globalMe or dr_dt, dr_beginStep
!!  are stored in the data fortran module for the Driver unit, Driver_data.
!!  The "dr_" is meant to indicate that the variable belongs to the Driver Unit.
!!  all other normally named variables i, j, etc are local variables.
!!
!!
!!***

#define WRITE_TO_TECPLOT 1

subroutine Driver_evolveFlash()

  use Driver_data, ONLY: dr_globalMe, dr_nbegin,                    &
                         dr_nend, dr_dt, dr_wallClockTimeLimit, &
                         dr_tmax, dr_simTime, dr_redshift,      &
                         dr_nstep, dr_dtOld, dr_dtNew,          &
                         dr_restart, dr_elapsedWCTime
  use Driver_interface, ONLY : Driver_sourceTerms, Driver_computeDt, &
       Driver_getElapsedWCTime
  use Logfile_interface,ONLY : Logfile_stamp, Logfile_close
  use Timers_interface, ONLY : Timers_start, Timers_stop, &
                               Timers_getSummary
  use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                                Grid_getListOfBlocks, &
                                Grid_getBlkIndexLimits
  use Grid_interface, only: Grid_getCellCoords, Grid_solvePoisson, GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN

  use IO_interface,      ONLY : IO_output,IO_outputFinal

  use IO_data , ONLY : IO_checkpointFileIntervalStep, io_plotFileNumber, IO_plotFileIntervalTime, IO_plotFileIntervalStep

  implicit none

#include "constants.h"
#include "Flash.h"
 include "Flash_mpi.h"

  integer :: localNumBlocks
  integer :: blockCount
  integer :: blockList(MAXBLOCKS)
  integer :: sweepDummy
  integer :: count, firstfileflag
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  logical :: endrun

  call Logfile_stamp( 'Entering evolution loop' , '[Driver_evolveFlash]')
  call Timers_start("evolution")

  count = 0
  dr_dt = 1
  dr_nstep = 0
  firstfileflag = 0
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  call Heat_AD(blockCount, blockList, dr_simTime, dr_dt, dr_dtOld, sweepDummy)

  count = 1
  dr_dt = 1 
  dr_nstep = 1
  firstfileflag = 0
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  call IO_output(dr_simTime,dr_dt,dr_nstep+1,dr_nbegin,endRun,PLOTFILE_ONLY)

  call Timers_stop("evolution")
  call Logfile_stamp( 'Exiting evolution loop' , '[Driver_evolveFlash]')
  call Timers_getSummary(dr_nstep)
  call Logfile_stamp( "FLASH run complete.", "LOGFILE_END")
  call Logfile_close()

  return
  
end subroutine Driver_evolveFlash



