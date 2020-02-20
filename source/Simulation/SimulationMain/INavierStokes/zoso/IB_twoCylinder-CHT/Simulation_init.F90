!!****if* source/Simulation/SimulationMain/INavierStokes/2D/IB_Cyl_parallelIBVP_HYPRE_VD_halfDiamBel/Simulation_init
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
                              sim_xMax, sim_yMax, sim_gCell, sim_waveA, sim_Tbulk, &
                              sim_sinkB, sim_invRe, sim_meshMe

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
 
  use IO_interface, ONLY :  IO_getScalar

  use Driver_data, only: dr_restart

  implicit none

  include 'Flash_mpi.h'
#include "constants.h"
#include "Flash.h"

  
  real :: xpt, ypt, dtheta, angle, tita, dsb
  integer :: i, b, ibd, nodelocpos, NumVertices, NumAelem

  real :: L1,L2,n1(MDIM),n2(MDIM)



  call RuntimeParameters_get('xmin',    sim_xMin)
  call RuntimeParameters_get('ymin',    sim_yMin)
  call RuntimeParameters_get('xmax',    sim_xMax)
  call RuntimeParameters_get('ymax',    sim_yMax)
  
  call RuntimeParameters_get('waveA',    sim_waveA)

  call RuntimeParameters_get('invRe',    sim_invRe)

end subroutine Simulation_init
