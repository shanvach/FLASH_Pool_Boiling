subroutine Simulation_init()
  
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

!===========================================================================

  call RuntimeParameters_get("smooth_radius", sim_smooth_radius)
  box_hwidth = int(sim_smooth_radius) + 0.5 

end subroutine Simulation_init

