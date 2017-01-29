!!****if* source/Simulation/SimulationMain/Flash2Convert/Driver_init
!!
!! NAME
!!  Driver_init
!!
!! SYNOPSIS
!!  Driver_init(integer (IN) :: myPE)
!!
!! DESCRIPTION
!!
!!  This version is specialized for the Flash 2 file converter setup.  
!!  Specifically this ensures that an old timestep (dtOld) is not set to prevent
!!  the file converter from crashing.
!!
!!  Perform the Driver unit initializations.
!!  Gets runtime parameters from the flash.par, or from checkpoint file
!!  if run is a restart.
!!  
!!  Initializes the simulation time, initial dt, timestep, whether 
!!  from scratch or restart.
!!  Also, verifies that the initial dt passes CFL criteria.
!!
!! ARGUMENTS
!!
!!  myPE - current processor
!!
!! PARAMETERS
!!  
!!   These are the runtime parameters used by the basic Driver unit.
!!   Your specific implementation may have more runtime parameters.
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory.
!!   You may have overwritten these values with the flash.par values
!!   for your specific run.  
!!
!!    dtinit [REAL]
!!        Initial timestep
!!    dtmax [REAL]
!!        Maximum timestep
!!    dtmin [REAL]
!!        Minimum timestep
!!    nbegin [INTEGER]
!!        First timestep
!!    nend [INTEGER]
!!        Maximum number of timesteps to take
!!    restart [BOOLEAN]
!!        Is this a restart run?
!!    tinitial [REAL]
!!        Initial simulation time
!!    tmax [REAL]
!!        Maximum simulation time
!!    wall_clock_time_limit [REAL]
!!        Total wall clock time limit (seconds)
!!    tstep_change_factor {REAL]
!!        factor to allow multiplicative increase in dt until it
!!        hits the CFL condition. This lets initial dt be
!!        very conservative initially, but increase rapidly to find the
!!        the optimum value.
!!
!!
!! NOTES
!!
!! The Driver unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. These variables begin with "dr_"
!! like, dr_globalMe or dr_dt, dr_beginStep, and are stored in FORTRAN
!! module Driver_data (in file Driver_data.F90). The other variables
!! are local to the specific routine and do not have the prefix "dr_"
!!
!!
!!***

subroutine Driver_init(myPE)
  use Driver_data
  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use IO_interface, ONLY : IO_getScalar  
  implicit none       

#include "Flash.h"
#include "constants.h" 
  
  integer, intent(in) :: myPE

  ! get the parameters needed by Driver
  call RuntimeParameters_get("nend", dr_nend)

  call RuntimeParameters_get("restart", dr_restart)
  call RuntimeParameters_get("tmax", dr_tmax)
  call RuntimeParameters_get("tinitial",dr_initialSimTime)

  call RuntimeParameters_get("dtinit",dr_dtInit)
  call RuntimeParameters_get("dtmin",dr_dtMin)
  call RuntimeParameters_get("dtmax",dr_dtMax)
  call RuntimeParameters_get("tstep_change_factor", dr_tstepChangeFactor)
  call RuntimeParameters_get("wall_clock_time_limit",dr_wallClockTimeLimit)
!!  This is not yet in use. This parameter controls timestep based
!!  upon the temperature.
!!  call RuntimeParameters_get("temp_factor", dr_tempFactor)
#ifdef GRAVITY
  !call RuntimeParameters_get("useGravity",dr_useGravity)
#else
  !dr_useGravity=.false.
#endif

  
  if (dr_restart) then
  
     !get values from the scalar list from the previous run
     call IO_getScalar("nstep", dr_nstep)
     call IO_getScalar("time", dr_simTime)
     call IO_getScalar("dt", dr_dt)
     !call IO_getScalar("dtOld", dr_dtOld)
     dr_dtOld = 0
     dr_nbegin = dr_nstep
     dr_initialSimTime = dr_simTime

  else if (.not. dr_restart) then
     call Driver_abortFlash("Flash2Converter must be run as a restart!")
  end if
  return

end subroutine Driver_init








