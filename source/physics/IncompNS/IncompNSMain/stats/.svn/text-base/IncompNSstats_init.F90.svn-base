

subroutine IncompNSstats_init()

  use Driver_data, only : dr_simTime  
  use IncompNS_data, only : ins_restart,ins_nstep
  use IncompNSstats_data, only : instats_IntervalStep,instats_restart_stats,&
                                 instats_n,instats_start_time
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

  ! How many steps between statistics:
  call RuntimeParameters_get("statistics_steps",instats_IntervalStep) 

  ! Should we restart the statistics?
  call RuntimeParameters_get("restart_statistics",instats_restart_stats)

  ! Compute initial time before start stats:
  call RuntimeParameters_get("start_statistics_time",instats_start_time)

  if (ins_restart .and. instats_restart_stats .and. (dr_simTime .ge. instats_start_time)) then

    ! Here estimate ensemble number. This overestimates the number if instats_start time > 0.
    instats_n = ins_nstep/instats_IntervalStep  

  else
    instats_n = 0
  endif 

  return

end subroutine

