!!****if* source/Simulation/SimulationMain/SolidMechanics/First_test_Rigid/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! ARGUMENTS
!!
!!   
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Initializes initial conditions for INS-isotropic turbulence problem.
!!
!!***

subroutine Simulation_init()

  use Simulation_data, ONLY : sim_xMin, sim_yMin, &
                              sim_xMax, sim_yMax, sim_gCell

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use sm_iointerface, only : sm_ioReadSolid
  use SolidMechanics_Data, only : sm_BodyInfo

  use ImBound_data, only : ib_nu

  implicit none

#include "constants.h"
#include "Flash.h"

  integer :: ibd,i

  call RuntimeParameters_get('xmin',    sim_xMin)
  call RuntimeParameters_get('ymin',    sim_yMin)
  call RuntimeParameters_get('xmax',    sim_xMax)
  call RuntimeParameters_get('ymax',    sim_yMax)

  sim_gCell = .true.

  ib_nu = 1.

end subroutine Simulation_init
