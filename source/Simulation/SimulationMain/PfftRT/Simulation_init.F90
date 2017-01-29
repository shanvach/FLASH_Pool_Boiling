subroutine Simulation_init()
  
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  
!===========================================================================

  call RuntimeParameters_get("smooth_step_min", sim_smooth_step_min)
  call RuntimeParameters_get("smooth_step_delta", sim_smooth_step_delta)
  call RuntimeParameters_get("smooth_step_max", sim_smooth_step_max)
  call RuntimeParameters_get("bufFact", sim_bufFact)

  call RuntimeParameters_get("writeflsmdata", sim_writeflsmdata)

  call RuntimeParameters_get("isolevel_1", isolevels(1))
  call RuntimeParameters_get("isolevel_2", isolevels(2))
  call RuntimeParameters_get("isolevel_3", isolevels(3))

  if ( sim_smooth_step_delta .eq. 0.0 ) then
    call Driver_abortFlash("Simulation_init: smooth_step_delta must .ne. 0.0")
  endif
  
end subroutine Simulation_init

