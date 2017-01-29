
subroutine Simulation_init(myPE)

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get


  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: myPE


  
  call RuntimeParameters_get('xmin',sim_xMin)
  call RuntimeParameters_get('ymin',sim_yMin)
  call RuntimeParameters_get('zmin',sim_zMin)
  call RuntimeParameters_get('xmax',sim_xMax)
  call RuntimeParameters_get('ymax',sim_yMax)
  call RuntimeParameters_get('zmax',sim_zMax)
  
  
end subroutine Simulation_init
