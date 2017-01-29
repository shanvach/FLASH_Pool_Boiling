!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testI/sim_setErrorBars
!!
!!  NAME 
!!
!!   sim_setErrorBars
!!
!!  SYNOPSIS
!!
!!   sim_setErrorBars ()
!!
!!  DESCRIPTION
!!
!!   This routine sets the error bars for the simulation. These were determined
!!   from previous runs.
!!
!!***

subroutine sim_setErrorBars ()

  use Simulation_data,             ONLY : sim_focalPercentErrorX,   &
                                          sim_focalPercentErrorXZ,  &
                                          sim_powerPercentErrorX,   &
                                          sim_powerPercentErrorXZ,  &
                                          sim_rayPexitPercentError, &
                                          sim_rayXexitPercentError, &
                                          sim_rayZexitPercentError, &
                                          sim_refinementLevel,      &
                                          sim_symmetryTolerance

  implicit none
!
!
!     ...Set the symmetry tolerance, which checks, if the results are expected
!        due to symmetry.
!
!
  sim_symmetryTolerance = 1.e-8
!
!
!     ...Set the upper bound of the focal point percentage errors. Since these
!        depend on the refinement level, each level has its own upper bound.
!        The ones with the extension X apply to those rays having the z-coordinate
!        equal to the tube's center z-coordinate. The ones with the extension XZ
!        apply to those rays which will feel the deflecting force in both x- and
!        z-directions.
!
!
  sim_focalPercentErrorX  (1) = 0.99
  sim_focalPercentErrorX  (2) = 3.28
  sim_focalPercentErrorX  (3) = 0.56
  sim_focalPercentErrorX  (4) = 0.21
  sim_focalPercentErrorX  (5) = 0.01
  sim_focalPercentErrorX  (6) = 0.06

  sim_focalPercentErrorXZ (1) = 7.04
  sim_focalPercentErrorXZ (2) = 4.56
  sim_focalPercentErrorXZ (3) = 0.77
  sim_focalPercentErrorXZ (4) = 1.12
  sim_focalPercentErrorXZ (5) = 0.99
  sim_focalPercentErrorXZ (6) = 0.41
!
!
!     ...The same but for the power deposition errors.
!
!
  sim_powerPercentErrorX  (1) = 0.99
  sim_powerPercentErrorX  (2) = 0.46
  sim_powerPercentErrorX  (3) = 0.11
  sim_powerPercentErrorX  (4) = 0.02
  sim_powerPercentErrorX  (5) = 0.01
  sim_powerPercentErrorX  (6) = 0.01

  sim_powerPercentErrorXZ (1) = 1.24
  sim_powerPercentErrorXZ (2) = 0.04
  sim_powerPercentErrorXZ (3) = 0.10
  sim_powerPercentErrorXZ (4) = 0.05
  sim_powerPercentErrorXZ (5) = 0.09
  sim_powerPercentErrorXZ (6) = 0.03
!
!
!     ...Convert the info from refinement level -> individual rays.
!
!
  sim_rayXexitPercentError (1) = sim_focalPercentErrorX  (sim_refinementLevel)     ! launched with z = 0
  sim_rayXexitPercentError (2) = sim_focalPercentErrorX  (sim_refinementLevel)     ! launched with z = 0
  sim_rayXexitPercentError (3) = 0.01                                              ! launched with x = 0
  sim_rayXexitPercentError (4) = 0.01                                              ! launched with x = 0
  sim_rayXexitPercentError (5) = sim_focalPercentErrorXZ (sim_refinementLevel)
  sim_rayXexitPercentError (6) = sim_focalPercentErrorXZ (sim_refinementLevel)
  sim_rayXexitPercentError (7) = sim_focalPercentErrorXZ (sim_refinementLevel)
  sim_rayXexitPercentError (8) = sim_focalPercentErrorXZ (sim_refinementLevel)

  sim_rayZexitPercentError (1) = 0.01
  sim_rayZexitPercentError (2) = 0.01
  sim_rayZexitPercentError (3) = sim_focalPercentErrorX  (sim_refinementLevel)
  sim_rayZexitPercentError (4) = sim_focalPercentErrorX  (sim_refinementLevel)
  sim_rayZexitPercentError (5) = sim_focalPercentErrorXZ (sim_refinementLevel)
  sim_rayZexitPercentError (6) = sim_focalPercentErrorXZ (sim_refinementLevel)
  sim_rayZexitPercentError (7) = sim_focalPercentErrorXZ (sim_refinementLevel)
  sim_rayZexitPercentError (8) = sim_focalPercentErrorXZ (sim_refinementLevel)

  sim_rayPexitPercentError (1) = sim_powerPercentErrorX  (sim_refinementLevel)
  sim_rayPexitPercentError (2) = sim_powerPercentErrorX  (sim_refinementLevel)
  sim_rayPexitPercentError (3) = sim_powerPercentErrorX  (sim_refinementLevel)
  sim_rayPexitPercentError (4) = sim_powerPercentErrorX  (sim_refinementLevel)
  sim_rayPexitPercentError (5) = sim_powerPercentErrorXZ (sim_refinementLevel)
  sim_rayPexitPercentError (6) = sim_powerPercentErrorXZ (sim_refinementLevel)
  sim_rayPexitPercentError (7) = sim_powerPercentErrorXZ (sim_refinementLevel)
  sim_rayPexitPercentError (8) = sim_powerPercentErrorXZ (sim_refinementLevel)
!
!
!     ...Ready!
!
!
  return
end subroutine sim_setErrorBars
