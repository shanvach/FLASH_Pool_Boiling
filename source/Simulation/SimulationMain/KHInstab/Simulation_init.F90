subroutine Simulation_init()
  use Simulation_data
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stampMessage

  implicit none

#include "constants.h"
#include "Flash.h"

  integer :: i
  integer :: geom


  call Driver_getMype(MESH_COMM, sim_meshMe)

  call RuntimeParameters_get('perturbDens', sim_perturbDens)
  call RuntimeParameters_get('perturbVelx', sim_perturbVelx)
  call RuntimeParameters_get('perturbTemp', sim_perturbTemp)
  call RuntimeParameters_get('ambientDens', sim_ambientDens)
  call RuntimeParameters_get('ambientVelx', sim_ambientVelx)
  call RuntimeParameters_get('ambientTemp', sim_ambientTemp)

  call RuntimeParameters_get('perturbRadius', sim_perturbRadius)

end subroutine Simulation_init
