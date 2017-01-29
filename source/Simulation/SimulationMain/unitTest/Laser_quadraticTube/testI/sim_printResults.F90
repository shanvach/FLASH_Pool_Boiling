!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testI/sim_printResults
!!
!!  NAME 
!!
!!   sim_printResults
!!
!!  SYNOPSIS
!!
!!   sim_printResults ()
!!
!!  DESCRIPTION
!!
!!   This routine prints the obtained ray exit coordinates and power (together with their
!!   percentage errors) to a .dat file. Only the master processor can make the printout, as it
!!   is the only processor having all ray exit data.
!!
!!***

subroutine sim_printResults ()

  use Simulation_data,             ONLY : sim_baseName,             &
                                          sim_globalMe,             &
                                          sim_nFocalRays,           &
                                          sim_powerDecayFactor,     &
                                          sim_rayPexit,             &
                                          sim_rayPexitPercentError, &
                                          sim_rayXexit,             &
                                          sim_rayXexitPercentError, &
                                          sim_rayZexit,             &
                                          sim_rayZexitPercentError, &
                                          sim_xw,                   &
                                          sim_zw

  implicit none

#include "Flash.h"
#include "constants.h"

  character (len = MAX_STRING_LENGTH) :: fileName

  integer  :: fileUnit
  integer  :: ray
  integer  :: ut_getFreeFileUnit

  real     :: rayP, rayPanalytic, rayPerror, rayPerrorBar
  real     :: rayX, rayXanalytic, rayXerror, rayXerrorBar
  real     :: rayZ, rayZanalytic, rayZerror, rayZerrorBar
!
!
!   ...Do the printout only on the master processor.
!
!
  if (sim_globalMe /= MASTER_PE) then
      return
  end if
!
!
!   ...Open the printout file.
!
!
  fileUnit = ut_getFreeFileUnit ()
  fileName = trim (sim_baseName) // "Results.dat"

  open (fileUnit, file = fileName)
!
!
!     ...Print out the rays exit x-coordinates.
!
!
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "                         RAY X EXIT (errors in %)"
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') " Ray         X exit         X exit analytic     X error    X error bar"
  write (fileUnit,'(a)') " ---------------------------------------------------------------------"

  do ray = 1, sim_nFocalRays

     rayX         = sim_rayXexit (ray)
     rayXanalytic = sim_xw
     rayXerror    = abs (rayX - sim_xw) * 100.0 / sim_xw
     rayXerrorBar = sim_rayXexitPercentError (ray)

     write (fileUnit,'(i3,2f20.12,2f12.5)') ray, rayX, rayXanalytic, rayXerror, rayXerrorBar

  end do
!
!
!     ...Print out the rays exit z-coordinates.
!
!
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "                         RAY Z EXIT (errors in %)"
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') " Ray         Z exit         Z exit analytic     Z error    Z error bar"
  write (fileUnit,'(a)') " ---------------------------------------------------------------------"

  do ray = 1, sim_nFocalRays

     rayZ         = sim_rayZexit (ray)
     rayZanalytic = sim_zw
     rayZerror    = abs (rayZ - sim_zw) * 100.0 / sim_zw
     rayZerrorBar = sim_rayZexitPercentError (ray)

     write (fileUnit,'(i3,2f20.12,2f12.5)') ray, rayZ, rayZanalytic, rayZerror, rayZerrorBar

  end do
!
!
!     ...Print out the rays exit power.
!
!
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "                         RAY P EXIT (errors in %)"
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') " Ray         P exit         P exit analytic     P error    P error bar"
  write (fileUnit,'(a)') " ---------------------------------------------------------------------"

  do ray = 1, sim_nFocalRays

     rayP         = sim_rayPexit (ray)
     rayPanalytic = sim_powerDecayFactor
     rayPerror    = abs (rayP - sim_powerDecayFactor) * 100.0 / sim_powerDecayFactor
     rayPerrorBar = sim_rayPexitPercentError (ray)

     write (fileUnit,'(i3,2f20.12,2f12.5)') ray, rayP, rayPanalytic, rayPerror, rayPerrorBar

  end do
!
!
!   ...Close the printout file.
!
!
  close (fileUnit)
!
!
!     ...Ready!
!
!
  return
end subroutine sim_printResults
