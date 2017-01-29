!!****if* source/Simulation/SimulationMain/magnetoHD/BierShock/Simulation_init
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
  use Logfile_interface, ONLY : Logfile_stamp
  use RadTrans_interface, ONLY : RadTrans_mgdSetBc
  
  implicit none

  call RuntimeParameters_get('sim_skewFactor' , sim_skewFactor )
  call RuntimeParameters_get('sim_rho' , sim_rho )
  call RuntimeParameters_get('sim_eint', sim_eint)
  call RuntimeParameters_get('sim_velx', sim_velx)
  call RuntimeParameters_get('sim_erad', sim_erad)

  call RuntimeParameters_get('eos_singleSpeciesZ', sim_singleSpeciesZ)
  call RuntimeParameters_get('ymin' , sim_ymin )  
  call RuntimeParameters_get('ymax' , sim_ymax )  

end subroutine Simulation_init
