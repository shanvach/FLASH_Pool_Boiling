
module Cool_data

!==============================================================================

  implicit none
#include "Flash.h"
! Data tables containing cooling curve(s)

  integer, save :: cl_meshMe, cl_numProcs 

  character(len=80), save     :: cl_fctn_file
  integer, parameter          :: cl_Nmax = 1000
  real, dimension(cl_Nmax), save :: cl_logT, cl_ne, cl_nH, &
       cl_nt, cl_logLnet, cl_logLnorm, &
       cl_logU, cl_logtau, cl_P12, cl_rho24, cl_Ci, cl_mubar, &
                               cl_dlogLdlogU, cl_dnedlogU, cl_dntdlogU
  integer, save               :: cl_N
  real, save, dimension(NSPECIES) :: cl_Xin, cl_Abar
  logical, save :: useCool

! Parameter passing for ODE integrator

  real, save :: cl_rho, cl_gamma, cl_smalle, cl_Boltzmann, cl_AMU

!==============================================================================

end module Cool_data
