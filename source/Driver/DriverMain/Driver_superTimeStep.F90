!!****if* source/Driver/DriverMain/Driver_superTimeStep
!!
!! NAME
!!
!!  Driver_superTimeStep
!!
!! SYNOPSIS
!!
!!  Driver_superTimeStep()
!!
!! DESCRIPTION
!!
!! This routine implements the super time steppping advancement algorithm
!! to overcome small diffusive time scales in explicit formulation.
!!
!! REFERENCES
!!
!! "Super-Time-Stepping Acceleration of Explicit Schemes for Parabolic Problems",
!! V. Alexiades, G. Amiez, and P. Gremaud, Com. Num. Meth. Eng, 1996
!! 
!!
!!***

subroutine Driver_superTimeStep(dt,nuSTS,nstepSTS,nstepTotalSTS,dt_subSTS)

  implicit none

#include "constants.h"
  !! Argument list -----------------------
  real, intent(IN)    :: dt,nuSTS
  integer, intent(IN) :: nstepSTS,nstepTotalSTS
  real, intent(OUT)   :: dt_subSTS
  !! -------------------------------------

  !! Calculate a substep dt_subSTS (tau_i)
  dt_subSTS = dt/((nuSTS - 1.0)*cos((2.*nstepSTS - 1.0)/nstepTotalSTS * 0.5*PI) + 1.0 + nuSTS)

  return
end subroutine Driver_superTimeStep
