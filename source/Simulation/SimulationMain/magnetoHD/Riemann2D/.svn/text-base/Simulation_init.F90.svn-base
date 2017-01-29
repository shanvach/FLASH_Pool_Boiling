!!****if* source/Simulation/SimulationMain/magnetoHD/Riemann2D/Simulation_init
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
!!  It calls RuntimeParameters_get routine for initialization.
!!  Initializes initial conditions for a 2D version of the BrioWu Riemann problem.
!!
!!***

subroutine Simulation_init()

  use Simulation_data

  use Driver_interface, ONLY : Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use Logfile_interface, ONLY : Logfile_stamp

  implicit none

#include "constants.h"
#include "Flash.h"

  

  call Driver_getMype(MESH_COMM, sim_meshMe)

  call RuntimeParameters_get('smallx', sim_smallx)
  call RuntimeParameters_get('smallp', sim_smallP)
  call RuntimeParameters_get('gamma',  sim_gamma)
  call RuntimeParameters_get('dens_1', sim_dens1)
  call RuntimeParameters_get('pres_1', sim_pres1)
  call RuntimeParameters_get('velx_1', sim_velx1)
  call RuntimeParameters_get('vely_1', sim_vely1)
  call RuntimeParameters_get('velz_1', sim_velz1)
  call RuntimeParameters_get('magx_1', sim_magx1)
  call RuntimeParameters_get('magy_1', sim_magy1)
  call RuntimeParameters_get('magz_1', sim_magz1)

  call RuntimeParameters_get('dens_2', sim_dens2)
  call RuntimeParameters_get('pres_2', sim_pres2)
  call RuntimeParameters_get('velx_2', sim_velx2)
  call RuntimeParameters_get('vely_2', sim_vely2)
  call RuntimeParameters_get('velz_2', sim_velz2)
  call RuntimeParameters_get('magx_2', sim_magx2)
  call RuntimeParameters_get('magy_2', sim_magy2)
  call RuntimeParameters_get('magz_2', sim_magz2)

  call RuntimeParameters_get('dens_3', sim_dens3)
  call RuntimeParameters_get('pres_3', sim_pres3)
  call RuntimeParameters_get('velx_3', sim_velx3)
  call RuntimeParameters_get('vely_3', sim_vely3)
  call RuntimeParameters_get('velz_3', sim_velz3)
  call RuntimeParameters_get('magx_3', sim_magx3)
  call RuntimeParameters_get('magy_3', sim_magy3)
  call RuntimeParameters_get('magz_3', sim_magz3)

  call RuntimeParameters_get('dens_4', sim_dens4)
  call RuntimeParameters_get('pres_4', sim_pres4)
  call RuntimeParameters_get('velx_4', sim_velx4)
  call RuntimeParameters_get('vely_4', sim_vely4)
  call RuntimeParameters_get('velz_4', sim_velz4)
  call RuntimeParameters_get('magx_4', sim_magx4)
  call RuntimeParameters_get('magy_4', sim_magy4)
  call RuntimeParameters_get('magz_4', sim_magz4)

  call RuntimeParameters_get('killdivb',  sim_killdivb)
  if (NDIM == 1) then
     sim_killdivb = .false.
  endif

  sim_gcell = .true.

end subroutine Simulation_init
