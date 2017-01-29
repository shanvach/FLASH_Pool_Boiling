module Simulation_data

  implicit none
#include "Flash.h"
#include "constants.h"

  real, save :: domain_volume
  real, save :: sim_bufFact
  real, save :: sim_smooth_step_delta
  integer, save :: sim_smooth_step_max, sim_smooth_step_min
  real, dimension(NDIM), parameter :: gwidth_rel = (/ 0.05, 0.1, 0.15 /)
  real, dimension(MDIM), save :: gauss_ctr
  real, dimension(2, MDIM), save :: original_region_bb
  integer, dimension(MDIM), save :: gauss_ctr_idx
  real, parameter :: xrefine_rel = 0.4
  logical, save :: sim_writeflsmdata
  real, dimension(3), save :: isolevels

end module Simulation_data
