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

  sm_pk_velocity     =  -1.471e-05*time**6 + 0.0004412*time**5 - & 
                          0.004785*time**4 + 0.02112*time**3 - & 
                           0.01789*time**2 - 0.009742*time - 0.9982
  
  sm_pk_acceleration =  -8.824e-05*time**5 + 0.002206*time**4 - &
                          0.01914*time**3 + 0.06337*time**2 - &
                           0.03579*time - 0.009742


  return

end subroutine sm_pk_getVelocity

! Linear Slope 0.1
!
!  sm_pk_velocity     =  -1.471e-05*time**6 + 0.0004412*time**5 - & 
!                          0.004785*time**4 + 0.02112*time**3 - & 
!                           0.01789*time**2 - 0.009742*time - 0.9982
!  
!  sm_pk_acceleration =  -8.824e-05*time**5 + 0.002206*time**4 - &
!                          0.01914*time**3 + 0.06337*time**2 - &
!                           0.03579*time - 0.009742

! Linear Slope 0.2
!
!  sm_pk_velocity     =  -2.941e-05*time**6 + 0.0008824*time**5 - & 
!                          0.00957*time**4 + 0.04225*time**3 - & 
!                           0.03579*time**2 - 0.01948*time - 0.9964
!  
!  sm_pk_acceleration =  0.0001765*time**5 + 0.004412*time**4 - &
!                          0.03828*time**3 + 0.1267*time**2 - &
!                           0.07158*time - 0.01948

! Quadradic Slope 0.1
!
!  sm_pk_velocity     =  -7.095e-07*time**6 + 3.073e-05*time**5 - & 
!                          0.0005286*time**4 + 0.004584*time**3 - & 
!                           0.01092*time**2 + 0.007002*time - 1
!  
!  sm_pk_acceleration =  -4.257e-06*time**5 + 0.0001536*time**4 - &
!                          0.002114*time**3 + 0.01375*time**2 - &
!                           0.02183*time + 0.007002

! Quadradic Slope 0.3
!
!  sm_pk_velocity     =  2.778e-05*time**6 - 0.0005192*time**5 + & 
!                          0.002906*time**4 + 0.0004108*time**3 - & 
!                           0.01528*time**2 + 0.01273*time - 1
!  
!  sm_pk_acceleration =  0.0001667*time**5 - 0.002596*time**4 + &
!                          0.01162*time**3 + 0.001233*time**2 - &
!                           0.03055*time + 0.01273

! Polynomial for Experimental 100 Degree Data
!
!  sm_pk_velocity = -3.526e-13*time**18 + 4.865e-11*time**17 - 3.083e-09*time**16 &
!                   + 1.188e-07*time**15 - 3.112e-06*time**14 + 5.86e-05*time**13 &
!                   - 0.0008183*time**12 + 0.00862  *time**11 - 0.06894 *time**10 &
!                   + 0.4182   *time**9  - 1.906    *time**8  + 6.419   *time**7  &
!                   - 15.51    *time**6  + 25.76    *time**5  - 27.55   *time**4 &
!                   + 17.25    *time**3  - 5.411    *time**2  + 0.6141  *time - 1.002
!
!  sm_pk_acceleration = -6.346e-12 *time**17 + 8.271e-10 *time**16 - 4.932e-08*time**15 &
!                       + 1.782e-06*time**14 - 4.357e-05 *time**13 + 0.0007617*time**12 &
!                       - 0.00982  *time**11 + 0.09482   *time**10 - 0.6894   *time**9  &
!                       + 3.764    *time**8  - 15.25     *time**7  + 44.93    *time**6  &
!                       - 93.07    *time**5  + 128.8     *time**4  - 110.2    *time**3  &
!                       + 51.74    *time**2  - 10.82     *time     + 0.6141

! Velocity 08
!
!  sm_pk_velocity = -1.536e-07*time**9 + 8.022e-06*time**8 - 0.0001656*time**7 &
!                  + 0.001664 *time**6 - 0.007562 *time**5 + 0.003236 *time**4 &
!                  + 0.08781  *time**3 - 0.2156   *time**2 + 0.1374   *time - 1.012
!
!  sm_pk_acceleration = -1.382e-06*time**8 + 6.418e-05*time**7 - 0.001159*time**6 &
!                      + 0.009987 *time**5 - 0.03781  *time**4 + 0.01295 *time**3 &
!                      + 0.2634   *time**2 - 0.4313   *time    + 0.1374

! Velocity09
!
!  sm_pk_velocity = 6.787e-05*time**5 - 0.002118*time**4 + 0.02196*time**3 &
!                  - 0.06933 *time**2 + 0.05757 *time**2 - 1
!
!  sm_pk_acceleration = 0.0003394*time**4 - 0.008474*time**3 + 0.06587*time**2 &
!                      - 0.1387  *time + 0.05757

! Velocity10
!
!  sm_pk_velocity = 1.831e-05*time**5 - 0.0006848*time**4 + 0.008622*time**3 &
!                  - 0.03501 *time**2 + 0.04064  *time**2 - 1.003
!
!  sm_pk_acceleration = 9.154e-05*time**4 - 0.002739*time**3 + 0.02587*time**2 &
!                      - 0.07002*time + 0.04064

! Velocity11
!
!  sm_pk_velocity = -1.312e-06*time**7 + 6.55e-05*time**6 - 0.001226*time**5 &
!                   + 0.01043 *time**4 - 0.03866 *time**3 + 0.05574 *time**2 &
!                   - 0.02023 *time    - 1.002
!
!  sm_pk_acceleration = -9.184e-06*time**6 + 0.000393*time**5 - 0.006129*time**4 &
!                       + 0.04172 *time**3 - 0.116   *time**2 + 0.1115  *time - 0.02023

! Constant Vel
!  const_vel = -10.0
!
!  sm_pk_velocity = const_vel !* ( -1.261e-17*time - 1 )
!
!  sm_pk_acceleration = 0.0 !-1.261e-17*time






