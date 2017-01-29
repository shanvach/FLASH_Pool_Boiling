!!****if* source/Simulation/SimulationMain/INavierStokes/2D/LidDrivenCavity/Simulation_init
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

  implicit none

#include "constants.h"
#include "Flash.h"

  integer :: ibd,i

  call RuntimeParameters_get('xmin',    sim_xMin)
  call RuntimeParameters_get('ymin',    sim_yMin)
  call RuntimeParameters_get('xmax',    sim_xMax)
  call RuntimeParameters_get('ymax',    sim_yMax)

  sim_gCell = .true.

!!$  ibd = 1
!!$
!!$  allocate(sm_BodyInfo(ibd))
!!$  call sm_ioReadSolid(ibd)
!!$
!!$  do i=1,sm_BodyInfo(ibd)%nnp
!!$
!!$     write(*,*) sm_BodyInfo(ibd)%x(i),sm_BodyInfo(ibd)%y(i),sm_BodyInfo(ibd)%z(i)
!!$
!!$  enddo

end subroutine Simulation_init
