!!****if* source/Simulation/SimulationMain/magnetoHD/LULI/Simulation_init
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
!!
!! DESCRIPTION
!!
!!  Initializes all the parameters needed for a particular simulation
!!
!!
!! ARGUMENTS
!!
!!  
!!
!! PARAMETERS
!!
!!***

subroutine Simulation_init()
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
                                          RuntimeParameters_mapStrToInt
  use PhysicalConstants_interface, ONLY  : PhysicalConstants_get
  use Logfile_interface, ONLY : Logfile_stamp
  use Driver_interface, ONLY : Driver_abortFlash, &
                               Driver_getComm
  
  implicit none

#include "constants.h"
#include "Flash.h"

  character(len=MAX_STRING_LENGTH) :: str

  call RuntimeParameters_get('geometry', str)
  call RuntimeParameters_mapStrToInt(str, sim_geometry)

  call RuntimeParameters_get('sim_targetRadius', sim_targetRadius)
  call RuntimeParameters_get('sim_targetHeight', sim_targetHeight)
  call RuntimeParameters_get('sim_targetOffset', sim_targetOffset)
  call RuntimeParameters_get('sim_targetZOffset', sim_targetZOffset)
  call RuntimeParameters_get('sim_skewFactor', sim_skewFactor)
  
  call RuntimeParameters_get('sim_rhoTarg', sim_rhoTarg)
  call RuntimeParameters_get('sim_teleTarg', sim_teleTarg)
  call RuntimeParameters_get('sim_tionTarg', sim_tionTarg)
  call RuntimeParameters_get('sim_tradTarg', sim_tradTarg)
  
  call RuntimeParameters_get('sim_rhoCham', sim_rhoCham)
  call RuntimeParameters_get('sim_teleCham', sim_teleCham)
  call RuntimeParameters_get('sim_tionCham', sim_tionCham)
  call RuntimeParameters_get('sim_tradCham', sim_tradCham)

  call RuntimeParameters_get('smallX', sim_smallX)

  call RuntimeParameters_get('sim_pulseLength', sim_pulseLength)
  call RuntimeParameters_get('sim_laserEnergy', sim_laserEnergy)
  call RuntimeParameters_get('sim_targetGeom', sim_targetGeom)
  call RuntimeParameters_get('sim_driverType', sim_driverType)
  call RuntimeParameters_get('sim_computeBiermann', sim_computeBiermann)
  call RuntimeParameters_get('sim_ndiv', sim_ndiv)
  call RuntimeParameters_get('sim_useMesh', sim_useMesh)
  call RuntimeParameters_get('sim_meshGeom', sim_meshGeom)

  call PhysicalConstants_get('electron charge', sim_qele)
  call PhysicalConstants_get('speed of light', sim_speedlt)
  call PhysicalConstants_get('Avogadro', sim_avo)

  call Driver_getComm(MESH_COMM, sim_meshComm)  

  !! magnetic fields
  call RuntimeParameters_get('sim_Bx', sim_Bx)
  call RuntimeParameters_get('sim_By', sim_By)
  call RuntimeParameters_get('sim_Bz', sim_Bz)

#if defined(DIVB_VAR)
  call RuntimeParameters_get('killdivb',sim_killdivb)
#endif


end subroutine Simulation_init
