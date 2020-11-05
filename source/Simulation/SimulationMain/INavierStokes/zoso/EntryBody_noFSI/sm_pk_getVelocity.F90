!!!!
! Subroutine: sm_pk_getVelocity.F90
! Author: Elizabeth Gregorio
!
! Purpose: to calculate what the velocity and acceleration should be at a given
! point in time. Will be called in sm_pk_prescribed, sm_el01_mapParticles and
! ib_advect.
!          
!!!!
subroutine sm_pk_getVelocity(time,sm_pk_velocity,sm_pk_acceleration)

#include "constants.h"
  implicit none
  real, intent(in)    :: time
  real, intent(out)   :: sm_pk_velocity, sm_pk_acceleration


  sm_pk_velocity     =  -1.471e-05*time**6 + 0.0004412*time**5 - & 
                          0.004785*time**4 + 0.02112*time**3 - & 
                           0.01789*time**2 - 0.009742*time - 0.9982
  
  sm_pk_acceleration =  -8.824e-05*time**5 + 0.002206*time**4 - &
                          0.01914*time**3 + 0.06337*time**2 - &
                           0.03579*time - 0.009742

  return

end subroutine sm_pk_getVelocity
