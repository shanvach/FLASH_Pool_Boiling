!!****if* source/Simulation/SimulationMain/magnetoHD/BubbleShock/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_init()
!!
!! ARGUMENTS
!!
!!  none
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Initializes initial conditions for OrszagTang problem.
!!
!!***

subroutine Simulation_init()

  use Simulation_data, ONLY : sim_xMin, sim_yMin, sim_zMin,  &
                              sim_xMax,  sim_yMax, sim_zMax, sim_gCell, &
                              sim_killdivb, sim_lposn, sim_smallX,      &
                              sim_gam1,sim_gam2, sim_int1,sim_int2,     &
                              sim_temp1,sim_temp2,sim_gammac1,sim_gammac2, &
                              sim_densT,sim_presT,sim_velxT,sim_velyT,&
                              sim_velzT,sim_magxT,sim_magyT,sim_magzT,&
                              sim_densB,sim_presB,sim_velxB,sim_velyB,&
                              sim_velzB,sim_magxB,sim_magyB,sim_magzB,&
                              sim_bubbleRadius,sim_bubbleXCtr,sim_bubbleYCtr,sim_bubbleZCtr, & 
                              sim_bubbleDensity, &
                              sim_meshMe

  use Driver_interface, ONLY : Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Multispecies_interface, ONLY : Multispecies_getProperty

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "Multispecies.h"

  
  real, dimension(EOS_NUM) :: eosData
  real, dimension(NSPECIES) :: eosMassFr
  integer :: vecLen=1

  call Driver_getMype(MESH_COMM, sim_meshMe)

  call RuntimeParameters_get('xmin',    sim_xMin)
  call RuntimeParameters_get('ymin',    sim_yMin)
  call RuntimeParameters_get('zmin',    sim_zMin)
  call RuntimeParameters_get('xmax',    sim_xMax)
  call RuntimeParameters_get('ymax',    sim_yMax)
  call RuntimeParameters_get('zmax',    sim_zMax)
  call RuntimeParameters_get('killdivb',sim_killdivb)

  call RuntimeParameters_get('densT',sim_densT)
  call RuntimeParameters_get('presT',sim_presT)
  call RuntimeParameters_get('velxT',sim_velxT)
  call RuntimeParameters_get('velyT',sim_velyT)
  call RuntimeParameters_get('velzT',sim_velzT)
  call RuntimeParameters_get('magxT',sim_magxT)
  call RuntimeParameters_get('magyT',sim_magyT)
  call RuntimeParameters_get('magzT',sim_magzT)

  call RuntimeParameters_get('densB',sim_densB)
  call RuntimeParameters_get('presB',sim_presB)
  call RuntimeParameters_get('velxB',sim_velxB)
  call RuntimeParameters_get('velyB',sim_velyB)
  call RuntimeParameters_get('velzB',sim_velzB)
  call RuntimeParameters_get('magxB',sim_magxB)
  call RuntimeParameters_get('magyB',sim_magyB)
  call RuntimeParameters_get('magzB',sim_magzB)

  call RuntimeParameters_get('lposn',sim_lposn)
  call RuntimeParameters_get('bubbleRadius',sim_bubbleRadius)
  call RuntimeParameters_get('bubbleXCtr',sim_bubbleXCtr)
  call RuntimeParameters_get('bubbleYCtr',sim_bubbleYCtr)
  call RuntimeParameters_get('bubbleZCtr',sim_bubbleZCtr)
  call RuntimeParameters_get('bubbleDensity',sim_bubbleDensity)

  sim_gCell = .true.


  !! Get the gamma values for the two materials
  call Multispecies_getProperty(FLD1_SPEC,GAMMA,sim_gam1)
  call Multispecies_getProperty(FLD2_SPEC,GAMMA,sim_gam2)

! set a dummy index to stuff arrays to call the EOS
! set mass fractions for the first fluid 
!!$  eosMassFr(1) = 1.-sim_smallX
!!$  eosMassFr(2) = sim_smallX
!!$  eosdata(EOS_DENS)=sim_densB
!!$  eosdata(EOS_PRES)=sim_presB
!!$  call Eos(MODE_DENS_PRES,vecLen,eosData,eosMassFr)
!!$  sim_temp1=eosData(EOS_TEMP)
!!$  sim_gammac1=eosData(EOS_GAMC)
!!$  sim_int1=eosData(EOS_EINT)
!!$
!!$
!!$  eosMassFr(1) = sim_smallX
!!$  eosMassFr(2) = 1.-sim_smallX
!!$  eosdata(EOS_DENS)=sim_densT
!!$  eosdata(EOS_PRES)=sim_presT
!!$  call Eos(MODE_DENS_PRES,vecLen,eosData,eosMassFr)
!!$  sim_temp2=eosData(EOS_TEMP)
!!$  sim_gammac2=eosData(EOS_GAMC)
!!$  sim_int2=eosData(EOS_EINT)


end subroutine Simulation_init
