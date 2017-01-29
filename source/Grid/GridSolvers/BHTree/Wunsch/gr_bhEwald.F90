!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhEwald
!!
!! NAME
!!
!!  gr_bhEwald
!!
!!
!! SYNOPSIS
!!
!!   real field = gr_bhEwald(
!!                           real(in)    :: x,
!!                           real(in)    :: y,
!!                           real(in)    :: z
!!        )
!!
!! DESCRIPTION
!!
!!   Interpolates in the Ewald field and returns its value for the point x,y,z.
!!
!! ARGUMENTS
!!
!!  x   - x-coordinate of the point where the Ewald field is determined
!!  y   - y-coordinate of the point where the Ewald field is determined
!!  z   - z-coordinate of the point where the Ewald field is determined
!!
!! RESULT
!!
!!  Value of the Ewald field in a given point.
!!
!! NOTES
!!
!!***



real function gr_bhEwald(x, y, z)

  use gr_bhData, ONLY: gr_bhEwaldFieldNx, gr_bhEwaldFieldNy, gr_bhEwaldFieldNz &
      & , gr_bhTreeEwald, gr_bhEwaldXMax, gr_bhEwaldYMax, gr_bhEwaldZMax
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  real,intent(in)    :: x, y, z
  real :: dx, dy, dz, p, q, r
  integer :: i, j, k, ip1, jp1, kp1, ll

  dx = gr_bhEwaldXMax / (gr_bhEwaldFieldNx-1)
  i = floor(x/dx)
  if ((i .LT. 0) .or. (i .GE. gr_bhEwaldFieldNx)) then
    call Driver_abortFlash("FATAL: Ewald field i out of limits")
  endif
  dy = gr_bhEwaldYMax / (gr_bhEwaldFieldNy-1)
  j = floor(y/dy)
  if ((j .LT. 0) .or. (j .GE. gr_bhEwaldFieldNy)) then
    call Driver_abortFlash("FATAL: Ewald field j out of limits")
  endif
  dz = gr_bhEwaldZMax / (gr_bhEwaldFieldNz-1)
  k = floor(z/dz)
  if ((k .LT. 0) .or. (k .GE. gr_bhEwaldFieldNz)) then
    call Driver_abortFlash("FATAL: Ewald field k out of limits")
  endif
  
  if (i .EQ. gr_bhEwaldFieldNx-1) then
    ip1 = i
  else 
    ip1 = i+1
  endif
  if (j .EQ. gr_bhEwaldFieldNy-1) then
    jp1 = j
  else 
    jp1 = j+1
  endif
  if (k .EQ. gr_bhEwaldFieldNz-1) then
    kp1 = k
  else 
    kp1 = k+1
  endif

  p = (x - i*dx)/dx
  q = (y - j*dy)/dy
  r = (z - k*dz)/dz
  
  gr_bhEwald = ( &
    & (1.0 - p)*(1.0 - q)*(1.0 - r) * gr_bhTreeEwald(i  , j  , k  ) + &
    & (1.0 - p)*(1.0 - q)*(      r) * gr_bhTreeEwald(i  , j  , kp1) + &
    & (1.0 - p)*(      q)*(1.0 - r) * gr_bhTreeEwald(i  , jp1, k  ) + &
    & (1.0 - p)*(      q)*(      r) * gr_bhTreeEwald(i  , jp1, kp1) + &
    & (      p)*(1.0 - q)*(1.0 - r) * gr_bhTreeEwald(ip1, j  , k  ) + &
    & (      p)*(1.0 - q)*(      r) * gr_bhTreeEwald(ip1, j  , kp1) + &
    & (      p)*(      q)*(1.0 - r) * gr_bhTreeEwald(ip1, jp1, k  ) + &
    & (      p)*(      q)*(      r) * gr_bhTreeEwald(ip1, jp1, kp1)   &
    & )

  return
end

