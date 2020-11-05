!!!!
! Subroutine: sm_pk_prescribed.F90
! Author: Elizabeth Gregorio
!
! Purpose: to allow for feeding in velocity data in order to match the 
!          kinematics recorded in experiments.
!
!!!!
subroutine sm_pk_prescribed(time,maxrestparams,paramcoord,vc,vcd,vcdd)

#include "constants.h"
  use sm_pk_interface, only: sm_pk_getVelocity  

  implicit none
  integer, intent(in) :: maxrestparams
  real, intent(in)    :: time, paramcoord(maxrestparams)
  real, intent(out)   :: vc, vcd, vcdd

  real :: xo
  real :: sm_pk_velocity, sm_pk_acceleration

  !print*,"body vel = ",vel

  call sm_pk_getVelocity(time, sm_pk_velocity, sm_pk_acceleration)

  ! Parameters: Are s.t. x(t) = xo + vel*t 
  xo = paramcoord(1)

  vc  = xo + vel*time
  vcd =           vel
  vcdd=           0.0


  return

end subroutine sm_pk_prescribed
