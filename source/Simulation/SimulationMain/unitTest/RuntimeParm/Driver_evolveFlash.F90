!!****if* source/Simulation/SimulationMain/unitTest/RuntimeParm/Driver_evolveFlash
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
!! This is a very simple version of the Driver_evolveFlash routine,
!! that is meant to be used exclusively with RuntimeParameters Unit
!! testing. There is no time advancement involved here.
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

subroutine Driver_evolveFlash()

  use Driver_data, ONLY:   dr_nbegin,  dr_restart, dr_initialSimTime
  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : &
    RuntimeParameters_getNumInt, RuntimeParameters_getNumReal, &
    RuntimeParameters_getNumStr, RuntimeParameters_getNumLog, RuntimeParameters_getAll, &
    RuntimeParameters_setPrev, RuntimeParameters_getPrev

  implicit none

#include "constants.h"
#include "Flash.h"

  integer :: nIntParms, nRealParms, nStrParms, nLogParms
  logical :: restartTest
  real,allocatable, dimension(:) :: values
  logical,allocatable,dimension(:) :: changed
  character(len=MAX_STRING_LENGTH),allocatable, dimension(:) :: names
  integer :: i
  integer,parameter ::  iOut = 2

  !Testing the RuntimeParameter Unit
  print *, "initialSimTime = ", dr_initialSimTime
  call RuntimeParameters_getNumInt(nIntParms)
  call RuntimeParameters_getNumReal(nRealParms)
  call RuntimeParameters_getNumStr(nStrParms)
  call RuntimeParameters_getNumLog(nLogParms)
  
  print *, nIntParms, nRealParms, nStrParms, nLogParms      

  allocate(values(nRealParms))
  allocate(names(nRealParms))
  allocate(changed(nRealParms))
  call RuntimeParameters_getAll(nRealParms, names, values,changed)

  do i=1, nRealParms
     print *, names(i), values(i)
  end do

  open(iOut,file='unitTest_0000')

  if(dr_nbegin == 1) then
     call Driver_abortFlash("failed over ride flash.par test")
  else
     write(iOut,'("override of flash.par test passed")')
  end if

  dr_restart = .true.

  call RuntimeParameters_setPrev("restart", dr_restart)
  call RuntimeParameters_getPrev("restart", restartTest)
    
  if (dr_restart .EQV. restartTest) then
     write(iOut, '("setPrev and getPrev tests passed")')
     write(iOut, '("restart = ", L1)') dr_restart
     write(iOut, '("restartTest = ", L1)') restartTest
  else
     call Driver_abortFlash("setPrev and/or getPrev tests failed")
  end if
  deallocate(values)
  deallocate(names)
  deallocate(changed)

  write(iOut,'("all results conformed with expected values.")')
  close(iOut)

  return

end subroutine Driver_evolveFlash

