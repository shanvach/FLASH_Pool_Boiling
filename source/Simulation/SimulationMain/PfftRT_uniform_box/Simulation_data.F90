module Simulation_data

  implicit none
#include "Flash.h"

  real, save :: sim_smooth_radius, box_hwidth
  real, dimension(NDIM), parameter :: gwidth_rel = (/ 0.1, 0.2, 0.3 /)
  real, parameter :: xrefine_rel = 0.4

end module Simulation_data
