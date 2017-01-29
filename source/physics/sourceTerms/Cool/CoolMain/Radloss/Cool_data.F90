
Module Cool_data

#include "constants.h"
#include "Flash.h"

  integer, save :: cl_meshMe, cl_numProcs 
  integer, save :: cl_h1specindex, cl_elecSpecIndex
  real, save, dimension(NSPECIES) :: cl_Xin, cl_Abar
  real, save :: cl_tradmin, cl_tradmax, cl_dradmin, cl_dradmax
  logical,save :: useCool
  character(len=4) :: cl_speciesNameH1, cl_speciesNameElec

end Module Cool_data
