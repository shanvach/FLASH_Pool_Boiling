!!****if* source/Simulation/SimulationMain/INavierStokes/3D/IB_sphere_mcHYPRE_VD_wFS/Simulation_init
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
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Initializes initial conditions for INS-isotropic turbulence problem.
!!
!!***

subroutine Simulation_init()

  use Grid_data, only : gr_meshMe

  use Driver_data, ONLY: dr_simTime

  use Driver_interface, ONLY : Driver_abortFlash

  use Simulation_data, ONLY : sim_xMin, sim_yMin, &
                              sim_xMax, sim_yMax, sim_gCell

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use ImBound_data

  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies

  implicit none

#include "constants.h"
#include "Flash.h"

  
  real :: xpt, ypt, dtheta, angle, tita, dsb
  integer :: i, b, ibd, nodelocpos, NumVertices, NumAelem

  real :: L1,L2,n1(MDIM),n2(MDIM)



  call RuntimeParameters_get('xmin',    sim_xMin)
  call RuntimeParameters_get('ymin',    sim_yMin)
  call RuntimeParameters_get('xmax',    sim_xMax)
  call RuntimeParameters_get('ymax',    sim_yMax)

  sim_gCell = .true.

  ! Cylinder setup IB variables:
  Ro = 0.5 !0.5
  omg= 0.0 !2.0
  freq_nat = 0.*0.196 ! Natural shedding frequency at Re=185
  freq_t   = 1.0*freq_nat ! Translational motion frequency.
  ao = 0.2*2.*Ro ! Translational motion amplitude.


end subroutine Simulation_init
