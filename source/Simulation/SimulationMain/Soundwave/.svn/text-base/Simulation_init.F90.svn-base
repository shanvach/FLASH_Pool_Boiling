!!****if* source/Simulation/SimulationMain/Soundwave/Simulation_init
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
!!
!!
!!***

subroutine Simulation_init()
  
  use Simulation_data
  use Driver_interface, ONLY : Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp
  implicit none
#include "Flash.h"

  
  call Driver_getMype(MESH_COMM, sim_meshMe)

  call RuntimeParameters_get('orientation', sim_orientation)
  call RuntimeParameters_get('wavelength', sim_wavelength)
  call RuntimeParameters_get('rho_init', sim_rhoInit)
  call RuntimeParameters_get('perturb_amp', sim_perturbAmp)
  call RuntimeParameters_get('cs', sim_cs)
  call RuntimeParameters_get('smallp', sim_smallP)
  call RuntimeParameters_get('smallx', sim_smallX) 
  
  call RuntimeParameters_get('gamma', sim_gamma)
  

  call Logfile_stamp( "initializing Sod problem",  &
       "[Simulation_init]")
     

end subroutine Simulation_init





