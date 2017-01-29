!!****if* source/Simulation/SimulationMain/NucOToRT/sim_setupOutputGrid
!!
!! NAME
!!
!!  sim_setupOutputGrid
!!
!!
!! SYNOPSIS
!!
!!  call sim_setupOutputGrid()
!!
!!
!! DESCRIPTION
!!
!!  Sets up data for the output grid of the
!!  nucOutput-to-radTrans converter.
!!
!! ARGUMENTS
!!
!!
!!***

subroutine sim_setupOutputGrid()
  
  use sim_outputGridData

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp
  implicit none
#include "Flash.h"


  call Logfile_stamp( "setting up NucOToRT output grid",  &
       "[sim_setupOutputGrid]")

  allocate(outGridData(iOgB0:iOgE0, jOgB0:jOgE0, kOgB0:kOgE0, nvarsOg))
  outGridData = 0.0             !zero-fill
  allocate(outLineData(nvarsOg))

  allocate(outGridNumCount(iOgB0:iOgE0, jOgB0:jOgE0, kOgB0:kOgE0))
  outGridNumCount = 0             !zero-fill

  allocate(outGridMappedVol(iOgB0:iOgE0, jOgB0:jOgE0, kOgB0:kOgE0))
  outGridMappedVol = 0.0             !zero-fill

  

end subroutine sim_setupOutputGrid
