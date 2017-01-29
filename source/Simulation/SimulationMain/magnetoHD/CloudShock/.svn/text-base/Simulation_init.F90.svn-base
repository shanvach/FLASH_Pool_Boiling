!!****if* source/Simulation/SimulationMain/magnetoHD/CloudShock/Simulation_init
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
!!
!! Reference: Dai & Woodward, JCP, 142:331--369, 1998
!!
!!***

subroutine Simulation_init()

  use Simulation_data, ONLY : sim_gamma, sim_xMin, sim_yMin, sim_zMin,  &
                              sim_xMax,  sim_yMax, sim_zMax, sim_gCell, &
                              sim_smallP,sim_smallX,                    &
                              sim_killdivb, sim_lposn,                  &
                              sim_dR,sim_pR,sim_uR,sim_vR,sim_wR,sim_bxR,sim_byR,sim_bzR,&
                              sim_dL,sim_pL,sim_uL,sim_vL,sim_wL,sim_bxL,sim_byL,sim_bzL,&
                              sim_cloudRadius,sim_cloudXCtr,sim_cloudYCtr,sim_cloudZCtr, & 
                              sim_cloudDensity, sim_Type,               &
                              sim_meshMe

  use Driver_interface, ONLY : Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "constants.h"
#include "Flash.h"

  

  call Driver_getMype(MESH_COMM, sim_meshMe)

  call RuntimeParameters_get('gamma',   sim_gamma)
  call RuntimeParameters_get('xmin',    sim_xMin)
  call RuntimeParameters_get('ymin',    sim_yMin)
  call RuntimeParameters_get('zmin',    sim_zMin)
  call RuntimeParameters_get('xmax',    sim_xMax)
  call RuntimeParameters_get('ymax',    sim_yMax)
  call RuntimeParameters_get('zmax',    sim_zMax)
#if defined(DIVB_VAR)
  call RuntimeParameters_get('killdivb',sim_killdivb)
#endif
  call RuntimeParameters_get('smallp',  sim_smallP)
  call RuntimeParameters_get('smallx',  sim_smallX)

  call RuntimeParameters_get('d_left',sim_dL)
  call RuntimeParameters_get('p_left',sim_pL)
  call RuntimeParameters_get('u_left',sim_uL)
  call RuntimeParameters_get('v_left',sim_vL)
  call RuntimeParameters_get('w_left',sim_wL)
#if defined(MAGX_VAR) && defined(MAGY_VAR) && defined(MAGZ_VAR)
  call RuntimeParameters_get('bx_left',sim_bxL)
  call RuntimeParameters_get('by_left',sim_byL)
  call RuntimeParameters_get('bz_left',sim_bzL)
#endif

  call RuntimeParameters_get('d_right',sim_dR)
  call RuntimeParameters_get('p_right',sim_pR)
  call RuntimeParameters_get('u_right',sim_uR)
  call RuntimeParameters_get('v_right',sim_vR)
  call RuntimeParameters_get('w_right',sim_wR)
#if defined(MAGX_VAR) && defined(MAGY_VAR) && defined(MAGZ_VAR)
  call RuntimeParameters_get('bx_right',sim_bxR)
  call RuntimeParameters_get('by_right',sim_byR)
  call RuntimeParameters_get('bz_right',sim_bzR)
#endif

  call RuntimeParameters_get('lposn',sim_lposn)
  call RuntimeParameters_get('cloudRadius',sim_cloudRadius)
  call RuntimeParameters_get('cloudXCtr',sim_cloudXCtr)
  call RuntimeParameters_get('cloudYCtr',sim_cloudYCtr)
  call RuntimeParameters_get('cloudZCtr',sim_cloudZCtr)
  call RuntimeParameters_get('cloudDensity',sim_cloudDensity)
  call RuntimeParameters_get('simType',sim_Type)

  sim_gCell = .true.

end subroutine Simulation_init
