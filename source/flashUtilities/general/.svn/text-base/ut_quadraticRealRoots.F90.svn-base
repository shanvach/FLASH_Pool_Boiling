!!****if* source/flashUtilities/general/ut_quadraticRealRoots
!!
!! NAME
!!
!!  ut_quadraticRealRoots
!!
!! SYNOPSIS
!!
!!  call ut_quadraticRealRoots (real,    intent (in)  :: q1,
!!                              real,    intent (in)  :: q0,
!!                              logical, intent (in)  :: rescale,
!!                              integer, intent (out) :: nRoot, 
!!                              real,    intent (out) :: Root1,
!!                              real,    intent (out) :: Root2)
!!
!! DESCRIPTION
!!
!!  Calculates all real roots of the quadratic polynomial:
!!
!!                 x^2 + q1 * x + q0
!!
!!  An option for rescaling the coefficients is provided for cases of very large
!!  initial coefficients. The code does not decide automatically, if rescaling is
!!  needed. The explicit rescaling option must be provided by the user. If some
!!  applications are known to lead to extremely large quadratic coefficients, then
!!  it is always safe to enforce rescaling. The cost of the rescaling is an extra
!!  square root and two divisions.
!!
!!  The code also deals with the case in which the discriminant is very close
!!  to zero, but the individual terms are large and opposite. In this case the
!!  code automatically sets the discriminant equal to zero and thus induces a
!!  degenerate real root.
!!
!! ARGUMENTS
!!
!!  q1      : coefficient of x term
!!  q0      : independent coefficient
!!  rescale : rescaling indicator
!!  nRoot   : number of real roots found and returned
!!  Root1   : the 1st real root (if any)
!!  Root2   : the 2nd real root (if any)
!!
!! NOTES
!!
!!***

subroutine ut_quadraticRealRoots (q1, q0, rescale,   nRoot, Root1, Root2)
  
  implicit none

  real,    intent (in)  :: q1, q0
  logical, intent (in)  :: rescale
  integer, intent (out) :: nRoot
  real,    intent (out) :: Root1, Root2

  real    :: a0, a1
  real    :: k, x, y, z

  real, parameter :: accuracy = 1.d-12
  real, parameter :: half     = 5.d-1
  real, parameter :: one      = 1.d0
  real, parameter :: zero     = 0.d0
!
!
!     ...Handle special cases.
!
!
  if (q0 == zero .and. q1 == zero) then

      nRoot = 1
      Root1 = zero

  else if (q0 == zero) then

      nRoot = 2
      Root1 = - q1
      Root2 = zero

  else if (q1 == zero) then

      if (q0 < zero) then
          x = sqrt (-q0)
          nRoot = 2
          Root1 = x
          Root2 = - x
      else
          nRoot = 0
      end if

  else
!
!
!     ...The general case. Do rescaling (if requested).
!
!
      if (rescale) then

          x = abs (q1)
          y = sqrt (abs (q0))

          if (x > y) then
              k  = x
              a1 = sign (one , q1)
              a0 = q0 / (x * x)
          else
              k  = y
              a1 = q1 / y
              a0 = sign (one , q0)
          end if

      else
          a1 = q1
          a0 = q0
      end if
!
!
!     ...Determine the real roots of the quadratic. Note, that either a1 or a0 might
!        have become equal to zero due to underflow. But both cannot be zero.
!
!
      x = a1 * a1
      nRoot = 0

      if (a0 /= zero) then

          y = a0 + a0
          y = y + y
          z = max (abs (x) , abs (y))
          y = x - y
          z = abs (y / z)

          if (z < accuracy) then                    ! this catches the cases where the discriminant
              y = zero                              ! is considered equal to zero, but the a1^2 and the 4a0
          end if                                    ! are large and opposite in magnitude.

          if (y > zero) then

              y = sqrt (y)
              y = - half * (a1 + sign (one,a1) * y)
              nRoot = 2
              Root1 = y                             ! 1st root from x^2 + a1 * x + a0
              Root2 = a0 / y                        ! 2nd root from x^2 + a1 * x + a0

          else if (y == zero) then

              nRoot = 1
              Root1 = - half * a1                   ! only root from x^2 + a1 * x + a0

          end if

      else

          nRoot = 2
          Root1 = - a1                              ! nonzero root from x^2 + a1 * x
          Root2 = zero                              ! zero root from x^2 + a1 * x

      end if
!
!
!     ...Rescale the roots (if needed).
!
!
      if (rescale) then

          if (nRoot == 1) then
              Root1 = Root1 * k
          else if (nRoot == 2) then
              Root1 = Root1 * k
              Root2 = Root2 * k
          end if

      end if

  end if
!
!
!     ...Ready!
!
!
  return
end subroutine ut_quadraticRealRoots
