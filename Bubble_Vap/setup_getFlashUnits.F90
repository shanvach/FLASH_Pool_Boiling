!!****f* object/setup_flashUnits
!!
!! NAME
!!
!!  setup_getFlashUnits
!!
!!
!! SYNOPSIS
!!
!!
!!  setup_getFlashUnits(unit_names)
!!
!!  setup_getFlashUnits(character())
!!
!!
!! DESCRIPTION
!!
!!  Return a character array of size NUM_UNITS containing
!!  the names of all of the FLASH units used to assemble
!!  the current executable
!!
!!  The unit_names variable should be declared as
!!
!!    use flashUnits
!!
!!  
!!    character (len=MAX_STRING_LENGTH) :: flash_units(NUM_UNITS) 
!!
!!
!!  The length of each character string is set to MAX_STRING_LENGTH,
!!  which is defined in the automatically generated flash_defines.fh
!!
!!***

  subroutine setup_getFlashUnits(unit_names)

#include "constants.h"
    implicit none

    integer, PARAMETER :: NUM_UNITS = 28
    character (len=MAX_STRING_LENGTH) :: unit_names(NUM_UNITS)
    character (len=80), PARAMETER :: all_unit_names(NUM_UNITS) = (/&
    &"Driver/DriverMain/INSfracstep                                                ",&
    &"Grid/GridBoundaryConditions                                                  ",&
    &"Grid/GridMain/paramesh/paramesh4/Paramesh4dev/PM4_package/headers            ",&
    &"Grid/GridMain/paramesh/paramesh4/Paramesh4dev/PM4_package/mpi_source         ",&
    &"Grid/GridMain/paramesh/paramesh4/Paramesh4dev/PM4_package/source             ",&
    &"Grid/GridMain/paramesh/paramesh4/Paramesh4dev/PM4_package/utilities/multigrid",&
    &"Grid/GridMain/paramesh/paramesh4/Paramesh4dev/flash_avoid_orrery             ",&
    &"Grid/GridSolvers/HYPRE_KPD/paramesh                                          ",&
    &"Grid/GridSolvers/MultigridMC_VarDens_HYPRE/poisson                           ",&
    &"Grid/localAPI                                                                ",&
    &"IO/IOMain/hdf5/serial/PM                                                     ",&
    &"IO/localAPI                                                                  ",&
    &"PhysicalConstants/PhysicalConstantsMain                                      ",&
    &"RuntimeParameters/RuntimeParametersMain                                      ",&
    &"Simulation/SimulationMain/INavierStokes/2D/Tecplot2D                         ",&
    &"Simulation/SimulationMain/INavierStokes/paramesh_routines/vardens_MG         ",&
    &"Simulation/SimulationMain/INavierStokes/zoso/Bubble_Vap                      ",&
    &"flashUtilities/contiguousConversion                                          ",&
    &"flashUtilities/general                                                       ",&
    &"flashUtilities/interpolation/oneDim                                          ",&
    &"flashUtilities/nameValueLL                                                   ",&
    &"flashUtilities/system/memoryUsage/legacy                                     ",&
    &"monitors/Logfile/LogfileMain                                                 ",&
    &"monitors/Timers/TimersMain/MPINative                                         ",&
    &"physics/Heat_AD/Heat_ADMain/vardens                                          ",&
    &"physics/IncompNS/IncompNSMain/extras                                         ",&
    &"physics/IncompNS/IncompNSMain/vardensTHERM                                   ",&
    &"physics/Multiphase/MultiphaseMain/INStherm                                   "/)

!!  if MAX_STRING_LENGTH < 80, then only the first MAX_STRING_LENGTH
!!  chars get copied over. So there is no problem with overflow here

    unit_names = all_unit_names 

    return

  end subroutine setup_getFlashUnits

  subroutine setup_getNumFlashUnits(numUnits)

#include "constants.h"
    implicit none

    integer, intent(out) :: numUnits
    integer, PARAMETER :: NUM_UNITS = 28

    numUnits = NUM_UNITS

    return

  end subroutine setup_getNumFlashUnits

