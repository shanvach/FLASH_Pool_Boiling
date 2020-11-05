!!!!
! Subroutine: sm_pk_getVelocity.F90
! Author: Elizabeth Gregorio
!
! Purpose: to calculate what the velocity and acceleration should be at a given
! point in time. Will be called in sm_pk_prescribed, sm_el01_mapParticles and
! ib_advect.
!          
!!!!
subroutine sm_pk_getVelocity(time,velocity,acceleration)

#include "constants.h"
  implicit none
  real, intent(in)    :: time
  real, intent(out)   :: sm_pk_velocity, sm_pk_acceleration


  sm_pk_velocity     = -2.918e-06*time**7 + 0.0001144*time**6 - 0.001822*time**5 + &
                   0.01513*time**4 - 0.06986*time**3 + 0.1761*time**2 - 0.1189*time - 1
  
  sm_pk_acceleration = -2.042e-05*time**6 + 0.0006863*time**5 - 0.009109*time**4 + &
                   0.06053*time**3 - 0.2096*time**2 + 0.3522*time - 0.1189

  return

end subroutine sm_pk_getVelocity
