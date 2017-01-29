!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testI/sim_doAnalysis
!!
!!  NAME 
!!
!!   sim_doAnalysis
!!
!!  SYNOPSIS
!!
!!   sim_doAnalysis (logical (out) :: perfect)
!!
!!  DESCRIPTION
!!
!!   This routine performs the laser quadratic tube simulation.
!!
!! ARGUMENTS
!!
!!  perfect : the success indicator for the simulation
!!
!!***

subroutine sim_doAnalysis (perfect)

  use Simulation_data,             ONLY : sim_globalComm,           &
                                          sim_globalMe,             &
                                          sim_nFocalRays,           &
                                          sim_nRaysMax,             &
                                          sim_powerDecayFactor,     &
                                          sim_rayPexit,             &
                                          sim_rayPexitPercentError, &
                                          sim_rayXexit,             &
                                          sim_rayXexitPercentError, &
                                          sim_rayZexit,             &
                                          sim_rayZexitPercentError, &
                                          sim_symmetryTolerance,    &
                                          sim_xw,                   &
                                          sim_XZtypeRays,           &
                                          sim_zw

  implicit none

#include "EnergyDeposition.h"
#include "Flash.h"
#include "constants.h"

  include "Flash_mpi.h"

  logical, intent (out) :: perfect

  integer  :: error
  integer  :: ray

  real     :: focalPercentError
  real     :: focalPercentErrorBar
  real     :: powerPercentError
  real     :: powerPercentErrorBar
  real     :: symmetryError

  real     :: Pdiff (1:sim_nRaysMax)
  real     :: Xdiff (1:sim_nRaysMax)
  real     :: Zdiff (1:sim_nRaysMax)
!
!
!     ...Calculate the percentage error between the actual and analytical solution
!        for the rays X exit coordinates and the rays exit power. Compare to the
!        errors bars. Only the X coordinates and power values are checked for rays
!        number 1 and 5. The other Z coordinate and all other rays are checked in the
!        symmetry section. For 2D simulations only rays 1 and 2 are checked.
!
!        The 'perfect' indicator will be set initially to true as a default on all
!        processors. Only the master processor perfroms the testing and has thus the
!        possibility to change its status to false. The logical status of this indicator
!        will be broadcast from the master to all processors after the analysis has been
!        performed.
!
!
  perfect = .true.

  if (sim_globalMe == MASTER_PE) then

      do ray = 1,sim_nFocalRays
         Xdiff (ray) = abs (sim_rayXexit (ray) - sim_xw)
         Zdiff (ray) = abs (sim_rayZexit (ray) - sim_zw)
         Pdiff (ray) = abs (sim_rayPexit (ray) - sim_powerDecayFactor)
      end do

      focalPercentErrorBar  = sim_rayXexitPercentError (1)
      powerPercentErrorBar  = sim_rayPexitPercentError (1)
      focalPercentError     = Xdiff (1) * 100.0 / sim_xw
      powerPercentError     = Pdiff (1) * 100.0 / sim_powerDecayFactor

      perfect = perfect .and. (focalPercentError  <= focalPercentErrorBar)
      perfect = perfect .and. (powerPercentError  <= powerPercentErrorBar)

      symmetryError = abs (Xdiff (1) - Xdiff (2))
      perfect = perfect .and. (symmetryError <= sim_symmetryTolerance)
      symmetryError = abs (Pdiff (1) - Pdiff (2))
      perfect = perfect .and. (symmetryError <= sim_symmetryTolerance)

      if (sim_XZtypeRays) then

          do ray = 3,4
             symmetryError = abs (Xdiff (1) - Zdiff (ray))
             perfect = perfect .and. (symmetryError <= sim_symmetryTolerance)
             symmetryError = abs (Pdiff (1) - Pdiff (ray))
             perfect = perfect .and. (symmetryError <= sim_symmetryTolerance)
          end do

          focalPercentErrorBar = sim_rayXexitPercentError (5)
          powerPercentErrorBar = sim_rayPexitPercentError (5)
          focalPercentError    = Xdiff (5)  * 100.0 / sim_xw
          powerPercentError    = Pdiff (5)  * 100.0 / sim_powerDecayFactor

          perfect = perfect .and. (focalPercentError  <= focalPercentErrorBar)
          perfect = perfect .and. (powerPercentError  <= powerPercentErrorBar)

          do ray = 6,8
             symmetryError = abs (Xdiff (5) - Xdiff (ray))
             perfect = perfect .and. (symmetryError <= sim_symmetryTolerance)
             symmetryError = abs (Zdiff (5) - Zdiff (ray))
             perfect = perfect .and. (symmetryError <= sim_symmetryTolerance)
             symmetryError = abs (Xdiff (5) - Zdiff (ray))
             perfect = perfect .and. (symmetryError <= sim_symmetryTolerance)
             symmetryError = abs (Zdiff (5) - Xdiff (ray))
             perfect = perfect .and. (symmetryError <= sim_symmetryTolerance)
             symmetryError = abs (Pdiff (5) - Pdiff (ray))
             perfect = perfect .and. (symmetryError <= sim_symmetryTolerance)
          end do

      end if

  end if
!
!
!     ...Broadcast the 'perfect' indicator.
!
!
  call MPI_Bcast (perfect,        &
                  1,              &
                  MPI_LOGICAL,    &
                  MASTER_PE,      &
                  sim_globalComm, &
                  error           )
!
!
!     ...Ready!
!
!
  return
end subroutine sim_doAnalysis
