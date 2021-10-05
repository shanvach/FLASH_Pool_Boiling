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
  real :: const_vel

  sm_pk_velocity = 9.125e-32*time**11 - 9.858e-28*time**10 + 4.648e-24*time**9 &
                  - 1.256e-20*time**8 + 2.143e-17*time**7  - 2.401e-14*time**6 &
                  + 1.776e-11*time**5 - 8.461e-09*time**4  + 2.436e-06*time**3 &
                  - 0.0003672*time**2 + 0.02104  *time     - 0.9998

  sm_pk_acceleration = 1.004e-30*time**10 - 9.858e-27*time**9 + 4.184e-23*time**8 &
                     - 1.005e-19*time**7  + 1.5e-16  *time**6 - 1.44e-13 *time**5 &
                     + 8.878e-11*time**4  - 3.384e-08*time**3 + 7.308e-06*time**2 &
                     - 0.0007344*time     + 0.02104 

  return

end subroutine sm_pk_getVelocity

! hydrophillic with 100 degree contact angle -- longer to reach term vel
!  sm_pk_velocity = 9.125e-32*time**11 - 9.858e-28*time**10 + 4.648e-24*time**9 &
!                  - 1.256e-20*time**8 + 2.143e-17*time**7  - 2.401e-14*time**6 &
!                  + 1.776e-11*time**5 - 8.461e-09*time**4  + 2.436e-06*time**3 &
!                  - 0.0003672*time**2 + 0.02104  *time     - 0.9998
!
!  sm_pk_acceleration = 1.004e-30*time**10 - 9.858e-27*time**9 + 4.184e-23*time**8 &
!                     - 1.005e-19*time**7  + 1.5e-16  *time**6 - 1.44e-13 *time**5 &
!                     + 8.878e-11*time**4  - 3.384e-08*time**3 + 7.308e-06*time**2 &
!                     - 0.0007344*time     + 0.02104 

! hydrophillic with 100 degree contact angle -- shorter to reach term vel
!  sm_pk_velocity = 1.018e-31*time**11 - 1.111e-27*time**10 + 5.3e-24*time**9 &
!                 - 1.449e-20*time**8  + 2.505e-17*time**7 - 2.843e-14*time**6 &
!                 + 2.128e-11*time**5  - 1.024e-08*time**4 + 2.966e-06*time**3 &
!                 - 0.0004479*time**2  + 0.02565  *time    - 0.9997
!
!  sm_pk_acceleration = 1.12e-30*time**10 - 1.111e-26*time**9 + 4.77e-23*time**8 &
!                     - 1.159e-19*time**7 + 1.753e-16*time**6 - 1.706e-13*time**5 &
!                     + 1.064e-10*time**4 - 4.096e-08*time**3 + 8.897e-06*time**2 &
!                     - 0.0008958*time    + 0.02565


