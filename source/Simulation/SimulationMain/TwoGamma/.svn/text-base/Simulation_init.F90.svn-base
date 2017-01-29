!!****if* source/Simulation/SimulationMain/TwoGamma/Simulation_init
!!
!! NAME
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!  call Simulation_init( )
!!
!!
!! DESCRIPTION
!!  Initializes all the parameters needed for the TwoGamma Problem
!!
!! ARGUMENTS
!!
!!  
!!
!!***

subroutine Simulation_init()
  
  use Simulation_data, ONLY : sim_xmin, sim_xmax, sim_ymin, sim_ymax, &
                              sim_small,sim_p0,sim_rho1,sim_rho2,&
                              sim_gam1,sim_gam2,sim_int1,sim_int2,&
                              sim_cvelx,sim_gammac1, sim_gammac2,sim_xpert,&
                              sim_temp1,sim_temp2
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Multispecies_interface, ONLY : Multispecies_getProperty
  use Eos_interface, ONLY : Eos
  implicit none

#include "Flash.h"
#include "Eos.h"
#include "Multispecies.h"
#include "constants.h"

  
  real, dimension(EOS_NUM) :: eosData
  real, dimension(NSPECIES) :: eosMassFr
  integer :: vecLen=1

  !! These parameters are defined by other units
  call RuntimeParameters_get( "xmin", sim_xmin)
  call RuntimeParameters_get( "xmax", sim_xmax)
  call RuntimeParameters_get( "ymin", sim_ymin)
  call RuntimeParameters_get( "ymax", sim_ymax)
  call RuntimeParameters_get( "small",sim_small)

  !! These parameters are defined by the this simulation
  call RuntimeParameters_get( "sim_p0", sim_p0)
  call RuntimeParameters_get( "sim_rho1", sim_rho1)
  call RuntimeParameters_get( "sim_rho2", sim_rho2)
  call RuntimeParameters_get( "sim_cvelx", sim_cvelx)
  
  !! Get the gamma values for the two materials
  call Multispecies_getProperty(FLD1_SPEC,GAMMA,sim_gam1)
  call Multispecies_getProperty(FLD2_SPEC,GAMMA,sim_gam2)


! set a dummy index to stuff arrays to call the EOS
! set mass fractions for the first fluid 
  eosMassFr(1) = 1.0e0 - sim_small
  eosMassFr(2) = sim_small
  eosdata(EOS_DENS)=sim_rho1
  eosdata(EOS_PRES)=sim_p0
  call Eos(MODE_DENS_PRES,vecLen,eosData,eosMassFr)
  sim_temp1=eosData(EOS_TEMP)
  sim_gammac1=eosData(EOS_GAMC)
  sim_int1=eosData(EOS_EINT)


  eosMassFr(1) = sim_small
  eosMassFr(2) = 1.0e0 - sim_small
  eosdata(EOS_DENS)=sim_rho2
  eosdata(EOS_PRES)=sim_p0
  call Eos(MODE_DENS_PRES,vecLen,eosData,eosMassFr)
  sim_temp2=eosData(EOS_TEMP)
  sim_gammac2=eosData(EOS_GAMC)
  sim_int2=eosData(EOS_EINT)

  !!  set the location of the interface
  sim_xpert = (sim_xmax - sim_xmin)/2.0

  return
end subroutine Simulation_init
