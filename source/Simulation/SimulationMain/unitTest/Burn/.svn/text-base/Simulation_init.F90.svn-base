!!****if* source/Simulation/SimulationMain/unitTest/Burn/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!! SYNOPSIS
!!
!!  call Simulation_init()
!!
!! DESCRIPTION
!!
!!  The routine initializes the unit scope data needed by the simulation unit.
!!  This includes runtime parameters and quantities derived from them
!!
!! ARGUMENTS
!!
!!   
!!
!!***
subroutine Simulation_init()

  use Simulation_data, ONLY:  sim_smallx, sim_imin, sim_imax, sim_jmin, sim_jmax, &
          sim_kmin, sim_kmax, sim_tempMin, sim_tempMax, sim_rhoMin, sim_rhoMax, sim_compA, sim_compB
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
 
!!  use RuntimeParameters_interface, ONLY:  RuntimeParameters_get

  implicit none


  
            
  call RuntimeParameters_get("smallx", sim_smallx)
  call RuntimeParameters_get("xmin", sim_imin)
  call RuntimeParameters_get("xmax", sim_imax)
  call RuntimeParameters_get("ymin", sim_jmin)
  call RuntimeParameters_get("ymax", sim_jmax)
  call RuntimeParameters_get("zmin", sim_kmin)
  call RuntimeParameters_get("zmax", sim_kmax)
  call RuntimeParameters_get("tempMin", sim_tempMin)
  call RuntimeParameters_get("tempMax", sim_tempMax)
  call RuntimeParameters_get("rhoMin", sim_rhoMin)
  call RuntimeParameters_get("rhoMax", sim_rhoMax)

  
  return
end subroutine Simulation_init


