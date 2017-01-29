!!****if* source/Simulation/SimulationMain/AdvectMassScalars/sim_pulseShape3D
!!
!! NAME
!!
!!  sim_pulseShap3De()
!!
!!
!! SYNOPSIS
!!
!!      sim_pulseShape3D (x, y, z, fctn, phase, val)
!!      sim_pulseShape3D (real, real, real, integer, real, real)
!!
!!
!! DESCRIPTION
!!
!!      Returns the value of the pulse function in the range [0,1] at position x.
!!
!!
!! ARGUMENTS
!!
!!     x               position (real)
!!     y               position (real)
!!     z               position (real)
!!     fctn            integer switch determining pulse shape
!!                          1 = constant
!!                          2 = square wave
!!                          3 = Gaussian
!!                          4 = sinusoid
!!                          5 = triangle
!!     phase             phase difference between density and tracers (real)
!!     val             the return value (real)
!!
!! NOTES
!!
!!***

subroutine sim_pulseShape3D (x, y, z, fctn, phase, val)

  !==============================================================================
  use Simulation_data
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  real, intent(IN)    :: x, y,z,  phase
  integer, intent(IN) :: fctn
  real, intent(OUT)   :: val

  real    :: dxphase, dyphase, dzphase, xx, yy, zz, radius

  !------------------------------------------------------------------------------

  ! Dongwook - note that phase = 0.
  dxphase = phase * (sim_xmax - sim_xmin)
  dyphase = phase * (sim_ymax - sim_ymin)
  dzphase = phase * (sim_zmax - sim_zmin)

  if ( fctn == 1 ) then

     !------------------------------------------------------------------------------
     ! constant density field

     val = 1.e0

  else if ( fctn == 2 ) then

     !------------------------------------------------------------------------------
     ! Square wave pulse shape

     xx = x - sim_xcent + dxphase
     yy = y - sim_ycent + dyphase
     zz = z - sim_zcent + dzphase

!!$     if ( xx > sim_xmax ) xx = xx - (sim_xmax - sim_xmin)
!!$     if ( xx < sim_xmin ) xx = xx + (sim_xmax - sim_xmin)
!!$     if ( yy > sim_ymax ) yy = yy - (sim_ymax - sim_ymin)
!!$     if ( yy < sim_ymin ) yy = yy + (sim_ymax - sim_ymin)
!!$     if ( zz > sim_zmax ) zz = zz - (sim_zmax - sim_zmin)
!!$     if ( zz < sim_zmin ) zz = zz + (sim_zmax - sim_zmin)

     if (       abs(xx) <= sim_width &
          .and. abs(sim_multid*yy) <= sim_width &
          .and. abs(sim_multid*zz) <= sim_width ) then
        val = 1.e0
     else
        val = 0.e0
     end if

  else if ( fctn == 3 ) then

     !------------------------------------------------------------------------------
     ! Gaussian wave pulse shape

     xx = x - (sim_xcent + dxphase)
     yy = y - (sim_ycent + dyphase)
     zz = z - (sim_zcent + dzphase)

!!$     val = exp( -( xx**2 &
!!$                  +sim_multid*yy**2 &
!!$                  +sim_multid*zz**2 &
!!$                 )/sim_width**2 )

     radius = sqrt( xx**2 + yy**2 + zz**2)
     if (radius <= sim_width) then
        val = 1.0e0
     else
        val = 0.0e0
     endif

  else if ( fctn == 4 ) then

     !------------------------------------------------------------------------------
     ! sinusoidal pulse shape

     val = (1.e0-0.5e0*sim_multid) * (            sin(sim_twoPi*(x+dxphase)) &
                                     + sim_multid*sin(sim_twoPi*(y+dyphase)) &
                                     + sim_multid*sin(sim_twoPi*(z+dyphase)))

  else if ( fctn == 5 ) then

     !------------------------------------------------------------------------------
     ! triangular pulse shape

     xx = x + dxphase

     if ( xx > sim_xmax ) xx = xx - (sim_xmax - sim_xmin)
     if ( xx < sim_xmin ) xx = xx + (sim_xmax - sim_xmin)

     if ( xx <= 0.5e0 ) then 
        val = -1.e0 + 4.e0*xx
     else 
        val = +1.e0 - 4.e0*(xx - 0.5e0)
     end if 

  else

     call Driver_abortFlash('[sim_pulseShape] ERROR: invalid pulse shape')

  end if

  return
end subroutine sim_pulseShape3D
