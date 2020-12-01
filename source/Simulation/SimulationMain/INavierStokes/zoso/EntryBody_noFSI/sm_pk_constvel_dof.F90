!
!
!!!!
subroutine sm_pk_constvel_dof(time,maxrestparams,paramcoord,vc,vcd,vcdd)

#include "constants.h"
  implicit none
  integer, intent(in) :: maxrestparams
  real, intent(in)    :: time, paramcoord(maxrestparams)
  real, intent(out)   :: vc, vcd, vcdd
  
  real :: xo, vel
  real :: sm_pk_velocity, sm_pk_acceleration

  call sm_pk_getVelocity(time, sm_pk_velocity, sm_pk_acceleration)

  ! Parameters: Are s.t. x(t) = xo + vel*t 
  xo = paramcoord(1)
  vel= paramcoord(2) 

  !print*,"body vel = ",vel

  vc  = xo + (vel*time) * sm_pk_velocity 
  vcd =           vel   * sm_pk_velocity
  vcdd=           sm_pk_acceleration

  return

end subroutine sm_pk_constvel_dof
