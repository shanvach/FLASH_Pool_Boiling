!!****if* source/Simulation/SimulationMain/magnetoHD/unitTest/ConstBier/Simulation_init
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
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Logfile_interface, ONLY : Logfile_stamp

  implicit none
#include "Flash.h"
#include "constants.h"

  call RuntimeParameters_get('sim_ptot', sim_ptot)
  call RuntimeParameters_get('sim_nele1', sim_nele1)
  call RuntimeParameters_get('sim_nele2', sim_nele2)
  call RuntimeParameters_get('sim_pele1', sim_pele1)
  call RuntimeParameters_get('sim_pele2', sim_pele2)
  call RuntimeParameters_get('eos_singleSpeciesA', sim_singleSpeciesA)
  call RuntimeParameters_get('eos_singleSpeciesZ', sim_singleSpeciesZ)
  call RuntimeParameters_get('xmin', sim_xmin)
  call RuntimeParameters_get('xmax', sim_xmax)
  call RuntimeParameters_get('ymin', sim_ymin)
  call RuntimeParameters_get('ymax', sim_ymax)

  call PhysicalConstants_get('Avogadro', sim_avogadro)
  call PhysicalConstants_get('Boltzmann', sim_boltzmann)  
  call PhysicalConstants_get('electron charge', sim_qele)
  call PhysicalConstants_get('speed of light', sim_speedlt)
end subroutine Simulation_init
