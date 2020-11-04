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
  implicit none
  integer, intent(in) :: maxrestparams
  real, intent(in)    :: time, paramcoord(maxrestparams)
  real, intent(out)   :: vc, vcd, vcdd

  real :: xo

  xo = paramcoord(1)

  !print*,"body vel = ",vel

  vcd = -2.918e-06*time**7 + 0.0001144*time**6 - 0.001822*time**5 +
          0.01513*time**4 - 0.06986*time**3 + 0.1761*time**2 - 0.1189*time - 1
  ! Parameters: Are s.t. x(t) = xo + vel*t 
  vc  = xo + vcd*time
  vcdd= -2.042e-05*time**6 + 0.0006863*time**5 - 0.009109*time**4 +
          0.06053*time**3 - 0.2096*time**2 + 0.3522*time - 0.1189

  return

end subroutine sm_pk_prescribed
