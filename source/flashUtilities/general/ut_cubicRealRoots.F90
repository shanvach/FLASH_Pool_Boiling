!!****if* source/flashUtilities/general/ut_cubicRealRoots
!!
!! NAME
!!
!!  ut_cubicRealRoots
!!
!! SYNOPSIS
!!
!!  call ut_cubicRealRoots (real,    intent (in)  :: c2,
!!                          real,    intent (in)  :: c1,
!!                          real,    intent (in)  :: c0,
!!                          logical, intent (in)  :: firstRootOnly,
!!                          integer, intent (out) :: nRoot, 
!!                          real,    intent (out) :: Root1,
!!                          real,    intent (out) :: Root2,
!!                          real,    intent (out) :: Root3)
!!
!! DESCRIPTION
!!
!!  Calculates all real roots of the cubic polynomial:
!!
!!                 x^3 + c2 * x^2 + c1 * x + c0
!!
!!  The first real root (which always exists) is obtained using an optimized
!!  Newton-Raphson scheme, except for near triply degenerate roots, where
!!  the analytical formula is being used. The other remaining real roots (if any)
!!  are obtained through deflation into a quadratic.
!!
!!  An option for returning just the first real root is provided, thus bypassing
!!  the possible deflation and quadratic solving. This option is useful for example
!!  when solving quartic equations. In that case only the first root of the
!!  subsidary cubic is needed.
!!
!!  The cubic root solver can handle any size of cubic coefficients and there is
!!  no danger of overflow due to proper rescaling of the cubic polynomial.
!!
!! ARGUMENTS
!!
!!  c2            : coefficient of x^2 term
!!  c1            : coefficient of x term
!!  c0            : independent coefficient
!!  firstRootOnly : if true, only the first real root is returned
!!  nRoot         : number of real roots found and returned
!!  Root1         : the 1st real root (always exists)
!!  Root2         : the 2nd real root (if any)
!!  Root3         : the 3rd real root (if any)
!!
!! NOTES
!!
!!***

subroutine ut_cubicRealRoots (c2, c1, c0, firstRootOnly,   nRoot, Root1, Root2, Root3)
  
  implicit none

  real,    intent (in)  :: c2, c1, c0
  logical, intent (in)  :: firstRootOnly
  integer, intent (out) :: nRoot
  real,    intent (out) :: Root1, Root2, Root3

  logical :: converged
  logical :: moreRealRoots
  logical :: newton
  logical :: rescale

  integer :: n

  real    :: a0, a1, a2
  real    :: k, s, t, u, x, y, z
  real    :: tolerance

  real, parameter :: accuracy    = 1.d-12

  real, parameter :: zero        = 0.d0
  real, parameter :: one100th    = 1.d-2
  real, parameter :: one27th     = 1.d0 / 27.d0
  real, parameter :: two27th     = 2.d0 / 27.d0
  real, parameter :: eleven27th  = 11.d0 / 27.d0
  real, parameter :: sixteen27th = 16.d0 / 27.d0
  real, parameter :: one9th      = 1.d0 / 9.d0
  real, parameter :: third       = 1.d0 / 3.d0
  real, parameter :: half        = 5.d-1
  real, parameter :: one         = 1.d0
  real, parameter :: two         = 2.d0

  real, parameter :: pi          = 3.1415926535897932384d0

  real, parameter :: p1          = 1.09574d0           !
  real, parameter :: q1          = 3.23900d-1          ! Newton-Raphson coeffs for class 1 and 2
  real, parameter :: r1          = 3.23900d-1          !
  real, parameter :: s1          = 9.57439d-2          !

  real, parameter :: p3          = 1.14413d0           !
  real, parameter :: q3          = 2.75509d-1          ! Newton-Raphson coeffs for class 3
  real, parameter :: r3          = 4.45578d-1          !
  real, parameter :: s3          = 2.59342d-2          !

  real, parameter :: q4          = 7.71845d-1          ! Newton-Raphson coeffs for class 4
  real, parameter :: s4          = 2.28155d-1          !

  real, parameter :: p51         = 8.78558d-1          !
  real, parameter :: p52         = 1.92823d-1          !
  real, parameter :: p53         = 1.19748d0           !
  real, parameter :: p54         = 3.45219d-1          !
  real, parameter :: q51         = 5.71888d-1          !
  real, parameter :: q52         = 5.66324d-1          !
  real, parameter :: q53         = 2.83772d-1          ! Newton-Raphson coeffs for class 5 and 6
  real, parameter :: q54         = 4.01231d-1          !
  real, parameter :: r51         = 7.11154d-1          !
  real, parameter :: r52         = 5.05734d-1          !
  real, parameter :: r53         = 8.37476d-1          !
  real, parameter :: r54         = 2.07216d-1          !
  real, parameter :: s51         = 3.22313d-1          !
  real, parameter :: s52         = 2.64881d-1          !
  real, parameter :: s53         = 3.56228d-1          !
  real, parameter :: s54         = 4.45532d-3          !
!
!
!     ...Handle special cases.
!
!
  if (c0 == zero .and. c1 == zero .and. c2 == zero) then

      nRoot = 1
      Root1 = zero
      return

  else if (c0 == zero .and. c1 == zero) then

      nRoot = 2
      Root1 = - c2
      Root2 = zero
      return

  else if (c0 == zero) then

      rescale = .true.

      call ut_quadraticRealRoots (c2,c1,rescale,  nRoot,x,y)

      if (nRoot == 2) then
          nRoot = 3
          Root1 = x
          Root2 = y
          Root3 = zero
      else if (nRoot == 1) then
          nRoot = 2
          Root1 = x
          Root2 = zero
      else
          nRoot = 1
          Root1 = zero
      end if

      return

  end if
!
!
!     ...The general case. Rescale cubic polynomial, such that largest absolute coefficient
!        is (exactly!) equal to 1.
!
!
  x = abs (c2)
  y = sqrt (abs (c1))
  z = abs (c0) ** third

  u = max (x,y,z)

  if (u == x) then

      k  = one / x
      t  = k * k
      a2 = sign (one , c2)
      a1 = c1 * t
      a0 = c0 * k * t

  else if (u == y) then

      k  = one / y
      a2 = c2 * k
      a1 = sign (one , c1)
      a0 = c0 * k * k * k

  else

      k  = one / z
      a2 = c2 * k
      a1 = c1 * k * k
      a0 = sign (one , c0)

  end if

  k = one / k
!
!
!     ...Set the best Newton-Raphson root estimates for the cubic. The easiest and most
!        robust conditions are checked first. The most complicated ones are last and only
!        done when absolutely necessary.
!
!
  if (a0 == one) then

      x = - p1 + q1 * a1 - a2 * (r1 - s1 * a1)

      newton = .true.
      moreRealRoots = abs (a1 + one) < accuracy .and. abs (a2 + one) < accuracy

  else if (a0 == - one) then

      x = p1 - q1 * a1 - a2 * (r1 - s1 * a1)

      newton = .true.
      moreRealRoots = abs (a1 + one) < accuracy .and. abs (a2 - one) < accuracy

  else if (a1 == one) then

      if (a0 > zero) then
          x = a0 * (- q4 - s4 * a2)
      else
          x = a0 * (- q4 + s4 * a2)
      end if

      newton = .true.
      moreRealRoots = .false.

  else if (a1 == - one) then

      y = - two27th
      y = y * a2
      y = y * a2 - third
      y = y * a2

      if (a0 < y) then
          x = + p3 - q3 * a0 - a2 * (r3 + s3 * a0)       ! + guess
      else
          x = - p3 - q3 * a0 - a2 * (r3 - s3 * a0)       ! - guess
      end if

      newton = .true.

      z = - eleven27th * a2
      y = z + sixteen27th
      z = z - sixteen27th

      moreRealRoots = (a0 <= y + accuracy) .and. (a0 >= z - accuracy)

  else if (a2 == one) then

      y = a1 - third
      z = a0 - one27th

      if (abs (y) < third * accuracy) then
          y = zero
      end if

      if (abs (z) < one27th * accuracy) then
          z = zero
      end if

      if (abs (y) < one100th .and. abs (z) < one100th) then

          s = - third * y              ! use analyic formula
          t = half * (s + z)
          u = t * t - s * s * s

          if (u < zero) then

              z = sqrt (s)
              u = acos (t / (s * z))
              x = - (z + z) * cos (third * u) - third

          else

              u = sqrt (u)
              u = abs (t) + u
              u = u ** third
              u = - sign (u,t)
              if (u /= zero) then
                  x = u + s / u - third
              else
                  x = u - third
              end if

          end if

          newton = .false.

      else

          y = third * a1 - two27th

          if (a1 <= third) then

              if (a0 > y) then
                  x = - p51 - q51 * a0 + a1 * (r51 - s51 * a0)   ! - guess
              else
                  x = + p52 - q52 * a0 - a1 * (r52 + s52 * a0)   ! + guess
              end if

          else

              if (a0 > y) then
                  x = - p53 - q53 * a0 + a1 * (r53 - s53 * a0)   ! <-1/3 guess
              else
                  x = + p54 - q54 * a0 - a1 * (r54 + s54 * a0)   ! >-1/3 guess
              end if

          end if

          newton = .true.

      end if

      z = a1 + a1 + a1
      y = one27th * (two - z)
      z = one9th  * (z + z + a1 - two)

      moreRealRoots = (a1 < third + accuracy) .and. (a0 <= y + accuracy) .and. (a0 >= z - accuracy)

  else if (a2 == - one) then

      y = a1 - third
      z = a0 + one27th

      if (abs (y) < third * accuracy) then
          y = zero
      end if

      if (abs (z) < one27th * accuracy) then
          z = zero
      end if

      if (abs (y) < one100th .and. abs (z) < one100th) then

          s = - third * y              ! use analyic formula
          t = half * (- s + z)
          u = t * t - s * s * s

          if (u < zero) then

              z = sqrt (s)
              u = acos (t / (s * z))
              x = (z + z) * cos (third *(u - pi)) + third

          else

              u = sqrt (u)
              u = abs (t) + u
              u = u ** third
              u = - sign (u,t)
              if (u /= zero) then
                  x = u + s / u + third
              else
                  x = u + third
              end if

          end if

          newton = .false.

      else

          y = two27th - third * a1

          if (a1 <= third) then

              if (a0 < y) then
                  x = + p51 - q51 * a0 - a1 * (r51 + s51 * a0)   ! +1 guess
              else
                  x = - p52 - q52 * a0 + a1 * (r52 - s52 * a0)   ! -1 guess
              end if

          else

              if (a0 < y) then
                  x = + p53 - q53 * a0 - a1 * (r53 + s53 * a0)   ! >1/3 guess
              else
                  x = - p54 - q54 * a0 + a1 * (r54 - s54 * a0)   ! <1/3 guess
              end if
                  
          end if

          newton = .true.

      end if

      z = a1 + a1 + a1
      y = one27th * (z - two)
      z = one9th  * (two - z - z - a1)

      moreRealRoots = (a1 < third + accuracy) .and. (a0 >= y - accuracy) .and. (a0 <= z + accuracy)

  end if
!
!
!     ...Perform Newton-Raphson iterations (if needed).
!
!
  converged = .false.

  do while (newton .and. .not.converged)

     tolerance = abs (x) * accuracy

     z = x + a2
     y = x + z
     z = z * x + a1
     y = y * x + z
     z = z * x + a0     

     u = x
     x = x - z / y
     converged = abs (u - x) <= tolerance

  end do

  nRoot = 1
  Root1 = x * k

  if (firstRootOnly) return
!
!
!     ...Forward / backward deflate rescaled cubic (if needed) to check for other real roots.
!        Deflation must be performed using the actual root of the cubic, not the rescaled root.
!        Otherwise deflation errors will be enhanced when undoing the rescaling on the extra
!        roots.
!
!
  if (moreRealRoots) then

      if (abs (x) > third) then
          x = one / Root1
          a0 = - c0 * x
          a1 = (a0 - c1) * x
      else
          y  = c1
          a1 = c2 + Root1
          a0 = y + a1 * Root1
      end if

      rescale = .true.

      call ut_quadraticRealRoots (a1,a0,rescale,  n,x,y)

      if (n == 2) then
          nRoot = 3
          Root2 = x
          Root3 = y
      else if (n == 1) then
          nRoot = 2
          Root2 = x
      end if

  end if

  if (nRoot == 1) return
!
!
!     ...All roots (with possible duplicates) have been generated. Order the roots
!        according to the number scale (most positive first, most negative last)
!        and eliminate duplicates. Two roots are considered identical, if the
!        magnitude of their difference is less than the accuracy attainable on any
!        one of the roots.
!
!
  if (nRoot == 2) then

      if (Root2 > Root1) then
          x = Root1
          Root1 = Root2
          Root2 = x
      end if

      z = abs (Root1 - Root2)
      
      if (z <= abs (Root1 * accuracy)) then
          nRoot = 1
      end if

  else

      if (Root2 > Root1) then
          x = Root1
          Root1 = Root2
          Root2 = x
      end if

      if (Root3 > Root1) then
          x = Root1
          Root1 = Root3
          Root3 = x
      else if (Root3 > Root2) then
          x = Root2
          Root2 = Root3
          Root3 = x
      end if

      z = abs (Root1 - Root2)
      
      if (z <= abs (Root1 * accuracy)) then
          nRoot = 2
          Root2 = Root3
      end if

      z = abs (Root1 - Root2)
      
      if (z <= abs (Root1 * accuracy)) then
          nRoot = 1
      end if

  end if
!
!
!     ...Ready!
!
!
  return
end subroutine ut_cubicRealRoots
