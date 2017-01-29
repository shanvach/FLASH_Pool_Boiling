!!****if* source/physics/sourceTerms/Cool/CoolMain/SutherlandDopita/cool_integrate
!!
!! NAME
!!
!!  cool_integrate
!!
!! SYNOPSIS
!!
!!  cool_integrate(real, intent(INOUT)  :: e,
!!                 real, intent(IN)     :: t,
!!                 real, intent(IN)     :: dt)
!!
!! DESCRIPTION
!!
!!  This uses a fourth-order Runge-Kutta method to integrate the Sutherland-
!!  Dopita cooling function and compute the new internal energy. This was
!!  handled by a third-order method provided by SVODE in Flash 2. 
!!
!! ARGUMENTS
!!
!!   e : The initial/final energy
!!
!!   t : The current time to solve for
!!
!!   dt : The current change in timestep
!!
!!
!!
!!
!!***


subroutine cool_integrate (e,t,dt)

  implicit none
  
  real, intent(INOUT)  :: e
  real, intent(IN)     :: t
  real, intent(IN)     :: dt

  real :: dt2, dt6, dtn
  real :: e1
  real :: dedt,dedt1,dedt2

  ! Use a Runge-Kutta method

  dt2 = 0.5*dt
  dt6 = dt/6.
  dtn = dt + dt2
  
  call cool_deriv(e,dedt)

  e1 = e + dt2*dedt

  call cool_deriv(e1,dedt1)
  
  e1 = e + dt2*dedt1

  call cool_deriv(e1,dedt2)

  e1 = e + dt*dedt2
  dedt2 = dedt1 + dedt2

  call cool_deriv(e1,dedt1)

  e = e + dt6*(dedt + dedt1 + 2.*dedt2)
  
  return

end subroutine cool_integrate
