
module IncompNSstats_data

  implicit none

  logical, save :: instats_restart_stats ! Restart statistics flag.
  integer, save :: instats_IntervalStep  ! Timesteps between stats computation.
  integer, save :: instats_n             ! Ensemble number.
  real,    save :: instats_start_time    ! Initial time before start of stats.

end module
