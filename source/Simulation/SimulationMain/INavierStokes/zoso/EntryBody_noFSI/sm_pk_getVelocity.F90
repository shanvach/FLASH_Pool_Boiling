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

  sm_pk_velocity = -0.002428 *time**7 + 0.04247 *time**6 - 0.2775 *time**5 &
                   + 0.8175  *time**4 - 1.022   *time**3 + 0.427  *time**2 &
                   + 0.006078*time    - 1.009

  sm_pk_acceleration = -0.017*time**6 + 0.2548  *time**5 - 1.387  *time**4 &
                       + 3.27*time**3 - 3.065   *time**2 + 0.8539 *time    &
                       + 0.006078

  return

end subroutine sm_pk_getVelocity

! very long -1 velocity (>20ndt)
!  sm_pk_velocity = 1.972*10**(-14)*time**8 - 8.18*10**(-12)*time**7 &
!                   + 1.361*10**(-9)*time**6 - 1.149*10**(-7)*time**5 &
!                   + 5.091*10**(-6)*time**4 - 0.0001084*time**3 &
!                   + 0.0009645*time**2       - 0.002502*time**1 - 1.001
!
!  sm_pk_acceleration = 1.578*10**(-13)*time**7 - 5.726*10**(-11)*time**6 &
!                       + 8.165*10**(-9)*time**5 - 5.746*10**(-7)*time**4 &
!                       + 2.036*10**(-5)*time**3 - 0.0003253*time**2 &
!                       + 0.001929*time - 0.002502 

! medium long -1 velocity (=15ndt)
!  sm_pk_velocity = -3.031e-13*time**8 + 9.6e-11*time**7 - 1.254e-08*time**6 &
!                   + 8.699e-07*time**5 - 3.414e-05*time**4 + 0.0007392*time**3 &
!                   - 0.007545*time**2 + 0.02677*time - 0.9999 
!
!  sm_pk_acceleration = -2.425e-12*time**7 + 6.72e-10*time**6 - 7.526e-08*time**5 &
!                       + 4.349e-06*time**4 - 0.0001366*time**3 + 0.002218*time**2 &
!                       - 0.01509*time + 0.02677 

! crashed
!  sm_pk_velocity = -3.01e-05*time**14  + 0.001161*time**13  - 0.01985*time**12 &
!                   + 0.1982*time**11   - 1.279   *time**10  + 5.582  *time**9 &
!                   - 16.75 *time**8    + 34.49   *time**7   - 47.96  *time**6 &
!                   + 43.65 *time**5    - 24.75   *time**4   + 8.065  *time**3 &
!                   - 1.31  *time**2    + 0.07649 *time      - 0.9999
!
!  sm_pk_acceleration = -0.0004214*time**13  + 0.01509*time**12  - 0.2382*time**11 &
!                       + 2.181   *time**10  - 12.79  *time**9   + 50.24 *time**8  &
!                       - 134     *time**7   + 241.4  *time**6   - 287.8 *time**5  &
!                       + 218.3   *time**4   - 98.99  *time**3   + 24.19 *time**2  &
!                       - 2.62    *time + 0.07649

! adjusted
!  sm_pk_velocity = -0.0001534*time**11  + 0.004001*time**10  - 0.0428*time**9 &
!                   + 0.2371  *time**8   - 0.6852  *time**7   + 0.7378*time**6 &
!                   + 1.091   *time**5   - 4.059   *time**4   + 4.428 *time**3 &
!                   - 2.005   *time**2   + 0.31    *time      - 1.004
!
!  sm_pk_acceleration = -0.001688*time**10  + 0.04001*time**9 - 0.3852*time**8 &
!                       + 1.897  *time**7   - 4.796  *time**6 + 4.427 *time**5 &
!                       + 5.457  *time**4   - 16.24  *time**3 + 13.28 *time**2 &
!                       - 4.009  *time      + 0.31

! works!
!  sm_pk_velocity = -0.002428 *time**7 + 0.04247 *time**6 - 0.2775 *time**5 &
!                   + 0.8175  *time**4 - 1.022   *time**3 + 0.427  *time**2 &
!                   + 0.006078*time    - 1.009
!
!  sm_pk_acceleration = -0.017*time**6 + 0.2548  *time**5 - 1.387  *time**4 &
!                       + 3.27*time**3 - 3.065   *time**2 + 0.8539 *time    &
!                       + 0.006078

! super long
!  sm_pk_velocity = -1.914e-09*time**13  + 1.737e-07*time**12 - 6.992e-06*time**11 &
!                   + 0.0001643*time**10 - 0.002495 *time**9  + 0.0256   *time**8  &
!                   - 0.1801   *time**7  + 0.8628   *time**6  - 2.734    *time**5  &
!                   + 5.413    *time**4  - 6.052    *time**3  + 3.287    *time**2  &
!                   - 0.6457   *time     - 0.9858
!
!  sm_pk_acceleration = -2.488e-08*time**12 + 2.084e-06*time**11 - 7.691e-05*time**10 &
!                       + 0.001643*time**9  - 0.02245  *time**8  + 0.2048   *time**7  &
!                       - 1.261   *time**6  + 5.177    *time**5  - 13.67    *time**4  &
!                       + 21.65   *time**3  - 18.15    *time**2  + 6.573    *time     &
!                       - 0.6457






