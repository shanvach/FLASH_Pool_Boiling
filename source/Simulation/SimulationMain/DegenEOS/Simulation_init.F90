!!****if* source/Simulation/SimulationMain/DegenEOS/Simulation_init
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
!!  It calls RuntimeParameters_get routine for initialization.
!!  Initializes initial conditions for BrioWu problem.
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp
  use Eos_interface, ONLY : Eos
  use Driver_interface, ONLY : Driver_getMype
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  
  real, dimension(EOS_NUM) :: eosData
  real, dimension(NSPECIES) :: eosMassFr
  integer :: vecLen=1

  call RuntimeParameters_get('posn',   sim_posn)
  call RuntimeParameters_get('dens0',  sim_densH)
  call RuntimeParameters_get('pres0',  sim_pres0)
  call RuntimeParameters_get('Atwood', sim_Atwood)
  call RuntimeParameters_get('Mach',   sim_Mach)
  call Driver_getMype(MESH_COMM, sim_meshMe)
  if (NDIM == 1) then
     sim_killdivb = .false.
  endif
  sim_gcell = .true.

  eosMassFr = 1.

  !! density high
  eosdata(EOS_DENS) = sim_densH
  eosdata(EOS_PRES) = sim_pres0
  eosdata(EOS_TEMP) = 5.e7
  call Eos(MODE_DENS_PRES,vecLen,eosData,eosMassFr)
  sim_tempH = eosData(EOS_TEMP)
  sim_gamcH = eosData(EOS_GAMC)
  sim_eintH = eosData(EOS_EINT)


  !! density low
  eosdata(EOS_DENS) = sim_densH*(1.-sim_Atwood)/(1.+sim_Atwood)
  eosdata(EOS_PRES) = sim_pres0
  eosdata(EOS_TEMP) = 8e9
  call Eos(MODE_DENS_PRES,vecLen,eosData,eosMassFr)
  sim_tempL = eosData(EOS_TEMP)
  sim_gamcL = eosData(EOS_GAMC)
  sim_eintL = eosData(EOS_EINT)

end subroutine Simulation_init
