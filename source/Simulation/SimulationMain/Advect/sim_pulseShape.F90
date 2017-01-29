!!****if* source/Simulation/SimulationMain/Advect/sim_pulseShape
!!
!! NAME
!!
!!  sim_pulseShape()
!!
!!
!! SYNOPSIS
!!
!!    sim_pulseShape (real (IN) :: x, 
!!                    integer(IN) :: fctn)
!!
!!
!! DESCRIPTION
!!
!! Finds the  shape of the density pulse for the advection problem
!! and returns the value of the pulse function at (scaled) position x.
!!
!!
!! ARGUMENTS
!!
!!     x -       The distance of the given point from the
!!              pulse midplane, scaled to the pulse width.
!!     fctn -     integer switch determining pulse shape
!!               1 = square wave
!!               2 = Gaussian
!!
!! NOTES
!!
!!     Currently supports two shapes, square wave (fcn .eq. 1) and 
!!      Gaussian pulse (fcn .ne. 1)
!!     Has same name, but different argument list from AdvectMassScalars version
!!
!!***

real function sim_pulseShape (x, fctn)

!==============================================================================

  implicit none
  
  real,intent(IN) :: x
  integer,intent(IN)  :: fctn
  
  ! Square wave pulse shape
  if (fctn .eq. 1) then
     
     if (abs(x) .le. 1.) then
        sim_pulseShape = 1.
     else
        sim_pulseShape = 0.
     endif
     
     !------------------------------------------------------------------
     
     ! Gaussian wave pulse shape
     
  else
     
     sim_pulseShape = exp( -x**2. )
     
  endif
  
  
  return
end function sim_pulseShape
