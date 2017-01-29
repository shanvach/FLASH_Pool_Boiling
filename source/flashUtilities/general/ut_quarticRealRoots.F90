!!****if* source/flashUtilities/general/ut_quarticRealRoots
!!
!! NAME
!!
!!  ut_quarticRealRoots
!!
!! SYNOPSIS
!!
!!  call ut_quarticRealRoots (real,    intent (in)  :: q3,
!!                            real,    intent (in)  :: q2,
!!                            real,    intent (in)  :: q1,
!!                            real,    intent (in)  :: q0,
!!                            logical, intent (in)  :: printInfo,
!!                            integer, intent (out) :: nRoot, 
!!                            real,    intent (out) :: Root1,
!!                            real,    intent (out) :: Root2,
!!                            real,    intent (out) :: Root3,
!!                            real,    intent (out) :: Root4)
!!
!! DESCRIPTION
!!
!!  Calculates all real roots of the quartic polynomial:
!!
!!                 x^4 + q3 * x^3 + q2 * x^2 + q1 * x + q0
!!
!!  An option for printing a detailed info about the intermediate stages in solving
!!  the quartic is available. Since the code has not yet been extensively tested,
!!  this enables a detailed check in case something went wrong and the roots obtained
!!  are not proper.
!!
!!  The quartic root solver can handle any size of quartic coefficients and there is
!!  no danger of overflow, due to proper rescaling of the quartic polynomial.
!!
!!  The order of the real roots is according to the number scale (largest positive
!!  first, largest negative last).
!!
!! ARGUMENTS
!!
!!  q3        : coefficient of x^3 term
!!  q2        : coefficient of x^2 term
!!  q1        : coefficient of x term
!!  q0        : independent coefficient
!!  printInfo : if true, detailed info will be printed about intermediate stages
!!  nRoot     : number of real roots found and returned
!!  Root1     : the 1st real root (if any)
!!  Root2     : the 2nd real root (if any)
!!  Root3     : the 3rd real root (if any)
!!  Root4     : the 4th real root (if any)
!!
!! NOTES
!!
!!***

subroutine ut_quarticRealRoots (q3, q2, q1, q0, printInfo,   nRoot, Root1, Root2, Root3, Root4)
  
  use Driver_interface,  ONLY : Driver_abortFlash

  implicit none

  real,    intent (in)  :: q3, q2, q1, q0
  logical, intent (in)  :: printInfo
  integer, intent (out) :: nRoot
  real,    intent (out) :: Root1, Root2, Root3, Root4

  logical :: complexRoots1, complexRoots2
  logical :: deflate2Cubic, deflate2Quadratic
  logical :: firstRootOnly
  logical :: rescale

  integer :: i,m,n
  integer :: nLarge

  real    :: a0, a1, a2, a3
  real    :: d0, d1, d2
  real    :: k, s, t, u, w, x, y, z

  real    :: Root (1:4)

  real, parameter :: accuracy    = 1.d-12
  real, parameter :: isZero      = 1.d-16

  real, parameter :: zero        = 0.d0
  real, parameter :: one256th    = 3.90625d-3
  real, parameter :: one16th     = 6.25d-2
  real, parameter :: three8th    = 3.75d-1
  real, parameter :: fourth      = 2.5d-1
  real, parameter :: third       = 1.d0 / 3.d0
  real, parameter :: half        = 5.d-1
  real, parameter :: one         = 1.d0
!
!
!     ...Start.
!
!
  if (printInfo) then
      write (*,'(a)'        ) ' ------------------------------------------------'
      write (*,'(a,es24.16)') ' initial quartic q3    = ',q3
      write (*,'(a,es24.16)') ' initial quartic q2    = ',q2
      write (*,'(a,es24.16)') ' initial quartic q1    = ',q1
      write (*,'(a,es24.16)') ' initial quartic q0    = ',q0
      write (*,'(a)'        ) ' ------------------------------------------------'
  end if
!
!
!     ...Handle special cases. Since the cubic solver handles all its
!        special cases by itself, we need to check only for two cases:
!
!            1) independent term is zero -> solve cubic and include
!               the zero root
!
!            2) the biquadratic case.
!
!
  if (q0 == zero) then

      Root (1) = zero

      firstRootOnly = .false.

      call ut_cubicRealRoots (q3,q2,q1,firstRootOnly,  n,x,y,z)

      if (n == 1) then
          Root (2) = x
      else if (n == 2) then
          Root (2) = x
          Root (3) = y
      else
          Root (2) = x
          Root (3) = y
          Root (4) = z
      end if

      nRoot = 1 + n

  else if (q3 == zero .and. q1 == zero) then               ! solve biquadratic

      rescale = .true.

      call ut_quadraticRealRoots (q2,q0,rescale,  n,x,y)   ! none of the roots can be 0

      if (n == 2) then

          if (x > zero .and. y > zero) then

              x = sqrt (x)
              y = sqrt (y)
              nRoot = 4
              Root (1) = x
              Root (2) = y
              Root (3) = - x
              Root (4) = - y

          else if (x > zero) then

              x = sqrt (x)
              nRoot = 2
              Root (1) = x
              Root (2) = - x

          else if (y > zero) then

              y = sqrt (y)
              nRoot = 2
              Root (1) = y
              Root (2) = - y

          end if

      else if (n == 1) then

          if (x > zero) then

              x = sqrt (x)
              nRoot = 2
              Root (1) = x
              Root (2) = - x

          end if

      else

          nRoot = 0
          return

      end if

  else
!
!
!     ...The general case. Rescale quartic polynomial, such that largest absolute coefficient
!        is (exactly!) equal to 1.
!
!
      a3 = q3
      a2 = q2
      a1 = q1
      a0 = q0

      w = abs (a3)
      x = sqrt (abs (a2))
      y = abs (a1) ** third
      z = abs (a0) ** fourth

      u = max (w,x,y,z)

      if (u == w) then

          k  = one / w
          t  = k * k
          a3 = sign (one , a3)
          a2 = a2 * t
          a1 = a1 * k * t
          a0 = a0 * t * t

      else if (u == x) then

          k  = one / x
          t  = k * k
          a3 = a3 * k
          a2 = sign (one , a2)
          a1 = a1 * k * t
          a0 = a0 * t * t

      else if (u == y) then

          k  = one / y
          t  = k * k
          a3 = a3 * k
          a2 = a2 * t
          a1 = sign (one , a1)
          a0 = a0 * t * t

      else

          k  = one / z
          t  = k * k
          a3 = a3 * k
          a2 = a2 * t
          a1 = a1 * k * t
          a0 = sign (one , a0)

      end if

      k = one / k

      if (printInfo) then
          write (*,'(a,es24.16)') ' rescaling factor      = ',k
          write (*,'(a)'        ) ' ------------------------------------------------'
          write (*,'(a,es24.16)') ' rescaled quartic q3   = ',a3
          write (*,'(a,es24.16)') ' rescaled quartic q2   = ',a2
          write (*,'(a,es24.16)') ' rescaled quartic q1   = ',a1
          write (*,'(a,es24.16)') ' rescaled quartic q0   = ',a0
          write (*,'(a)'        ) ' ------------------------------------------------'
      end if
!
!
!     ...Form depressed quartic.
!
!
      if (a3 == zero) then

          w  = zero
          d2 = a2
          d1 = a1
          d0 = a0

      else if (a3 == one) then

          w = fourth
          x = three8th - a2
          y = one16th  - a1
          z = one256th - a0

          d2 = - x
          d1 = half * x - y
          y  = fourth * y - one16th * x
          d0 = y - z

      else if (a3 == - one) then

          w = - fourth
          x = three8th - a2
          y = one16th  + a1
          z = one256th - a0

          d2 = - x
          d1 = y - half * x
          y  = fourth * y - one16th * x
          d0 = y - z

      else

          w  = fourth * a3
          x  = w * w
          y  = x + x
          z  = x + y - a2
          t  = z + z

          d2 = - (t + a2)
          d1 = w * (t + y) + a1
          d0 = a0 - a1 * w - x * z

      end if

      if (abs (d2) < isZero) d2 = zero
      if (abs (d1) < isZero) d1 = zero
      if (abs (d0) < isZero) d0 = zero

      if (printInfo) then
          write (*,'(a,es24.16)') ' depressed quartic d2  = ',d2
          write (*,'(a,es24.16)') ' depressed quartic d1  = ',d1
          write (*,'(a,es24.16)') ' depressed quartic d0  = ',d0
          write (*,'(a)'        ) ' ------------------------------------------------'
      end if
!
!
!     ...If the independent term of the depressed quartic is zero,
!        find the roots by calling the cubic solver.
!
!
      complexRoots1 = .false.
      complexRoots2 = .false.

      if (d0 == zero) then

          Root (1) = zero

          firstRootOnly = .false.

          call ut_cubicRealRoots (zero,d2,d1,firstRootOnly,  n,x,y,z)

          if (n == 1) then
              Root (2) = x
          else if (n == 2) then
              Root (2) = x
              Root (3) = y
          else
              Root (2) = x
              Root (3) = y
              Root (4) = z
          end if

          nRoot = 1 + n
!
!
!     ...If necessary, enter the resolvent cubic section.
!
!
      else

          s = - d2
          t = - (d0 + d0 + d0 + d0)
          u = s * t - d1 * d1

          if (printInfo) then
              write (*,'(a,es24.16)') ' resolvent cubic s     = ',s
              write (*,'(a,es24.16)') ' resolvent cubic t     = ',t
              write (*,'(a,es24.16)') ' resolvent cubic u     = ',u
              write (*,'(a)'        ) ' ------------------------------------------------'
          end if

          firstRootOnly = .true.

          call ut_cubicRealRoots (s,t,u,firstRootOnly,  n,x,y,z)

          if (printInfo) then
              write (*,'(a,es24.16)') ' Cubic resolvent root  = ',x
              write (*,'(a)'        ) ' ------------------------------------------------'
          end if
!
!
!     ...Proceed with the quadratic decomposition of the depressed quartic:
!
!            x^4 + d2 * x^2 + d1 * x + d0  =  (x^2 + s * x + t) * (x^2 - s * x + u)
!
!
          s = fourth * x * x
          u = max (abs (s) , abs (d0))
          s = s - d0
          y = abs (s / u)

          u = max (abs (x) , abs (d2))
          t = x - d2
          z = abs (t / u)

          if (y > z) then

              s = sqrt (s)

              if (x < zero) then
                  t = half * x - s
                  u = d0 / t
              else
                  u = half * x + s
                  t = d0 / u
              end if

              s = d1 / (s + s)

          else

              s = sqrt (t)

              if (x * d1 < zero) then
                  t = half * x - d1 / (s + s)
                  u = d0 / t
              else
                  u = half * x + d1 / (s + s)
                  t = d0 / u
              end if

          end if

          if (printInfo) then
              write (*,'(a,es24.16)') ' quadratic factor s    = ',s
              write (*,'(a,es24.16)') ' quadratic factor t    = ',t
              write (*,'(a,es24.16)') ' quadratic factor u    = ',u
              write (*,'(a)'        ) ' ------------------------------------------------'
              write (*,'(a,es24.16)') ' recalculated d2       = ',- s * s + t + u
              write (*,'(a,es24.16)') ' recalculated d1       = ',- s * t + s * u
              write (*,'(a,es24.16)') ' recalculated d0       = ',t * u
              write (*,'(a)'        ) ' ------------------------------------------------'
          end if
!
!
!     ...Calculate all real roots of the two quadratic factors:
!
!            (x^2 + s * x + t)   and   (x^2 - s * x + u)
!
!
          rescale = .false.

          call ut_quadraticRealRoots (s,t,rescale,  n,x,y)

          if (n == 2) then
              Root (1) = x
              Root (2) = y
          else if (n == 1) then
              Root (1) = x
          else
              complexRoots1 = .true.
          end if

          call ut_quadraticRealRoots (-s,u,rescale,  m,x,y)

          if (m == 2) then
              Root (n + 1) = x
              Root (n + 2) = y
          else if (m == 1) then
              Root (n + 1) = x
          else
              complexRoots2 = .true.
          end if

          nRoot = n + m

      end if
!
!
!     ...Apply the shift on all roots and check, if any of the roots magnitude is
!        large (>= 1/4). If yes, we know that this rescaled root is accurate. Find the largest
!        of these roots and unscale that root for backward deflation of the original quartic polynomial
!        into a cubic. If none of the roots is >= 1/4, then we will deflate to a quadratic
!        using the complex quadratic factor with the largest magnitude roots.
!
!
      Root (1:nRoot) = Root (1:nRoot) - w

      if (printInfo) then
          do n = 1,nRoot
             write (*,'(a,es24.16)') ' Rescaled quartic root = ',Root (n)
          end do
          write (*,'(a)') ' ------------------------------------------------'
      end if

      z = zero

      nLarge = 0

      do n = 1,nRoot
         x = Root (n)
         y = abs (x)
         if (y >= fourth) then
             nLarge = nLarge + 1
             if (y > z) then
                 Root1 = x
                 z = y
             end if
         end if
      end do

      deflate2Cubic     = (nLarge >  0) .and. (nLarge < nRoot)
      deflate2Quadratic = (nLarge == 0) .and. (complexRoots1 .or. complexRoots2)
!
!
!     ...Deflate original quartic from highest magnitude root using backward deflation
!        (evaluate constant term first).
!
!
      if (deflate2Cubic) then

          Root1 = Root1 * k

          if (printInfo) then
              write (*,'(a,es24.16)') ' Largest quartic root  = ',Root1
              write (*,'(a)'        ) ' ------------------------------------------------'
          end if

          x = one / Root1
          u = - q0 * x
          t = (u - q1) * x
          s = (t - q2) * x

          if (printInfo) then
              write (*,'(a,es24.16)') ' Deflated cubic s      = ',s
              write (*,'(a,es24.16)') ' Deflated cubic t      = ',t
              write (*,'(a,es24.16)') ' Deflated cubic u      = ',u
              write (*,'(a)'        ) ' ------------------------------------------------'
          end if

          firstRootOnly = .false.

          call ut_cubicRealRoots (s,t,u,firstRootOnly,  n,x,y,z)

          Root (1) = Root1

          if (n == 1) then
              Root (2) = x
          else if (n == 2) then
              Root (2) = x
              Root (3) = y
          else
              Root (2) = x
              Root (3) = y
              Root (4) = z
          end if

          nRoot = 1 + n

      else if (deflate2Quadratic) then

          if (complexRoots1 .and. complexRoots2) then
              t = t + w * (w + s)
              u = u + w * (w - s)
              if (abs (t) > abs (u)) then
                  s = w + w + s
              else
                  t = u
                  s = w + w - s
              end if
          else if (complexRoots1) then
              t = t + w * (w + s)
              s = w + w + s
          else if (complexRoots2) then
              t = u + w * (w - s)
              s = w + w - s
          end if

          t = t * k * k                     ! bring quadratic coeffs from rescaled
          s = s * k                         ! form to original form

          if (printInfo) then
              write (*,'(a,es24.16)') ' Complex quadratic s   = ',s
              write (*,'(a,es24.16)') ' Complex quadratic t   = ',t
              write (*,'(a)'        ) ' ------------------------------------------------'
          end if

          x = one / t
          t = q0 * x
          s = (q1 - s * t) * x

          if (printInfo) then
              write (*,'(a,es24.16)') ' Deflated quadratic s  = ',s
              write (*,'(a,es24.16)') ' Deflated quadratic t  = ',t
              write (*,'(a)'        ) ' ------------------------------------------------'
          end if

          rescale = .true.

          call ut_quadraticRealRoots (s,t,rescale,  n,x,y)

          if (n == 2) then
              Root (1) = x
              Root (2) = y
          else if (n == 1) then
              Root (1) = x
          end if

          nRoot = n

      end if

  end if
!
!
!     ...All roots (with possible duplicates) have been generated. Order the roots
!        according to the number scale (most positive first, most negative last)
!        and eliminate duplicates. Two roots are considered identical, if the
!        magnitude of their difference is less than the accuracy attainable on any
!        one of the roots.
!
!
  if (nRoot > 1) then

      do n = 1,nRoot-1                          ! ordering of roots
         i = n
         x = Root (n)
         do m = n+1,nRoot
            y = Root (m)
            if (y > x) then
                i = m
                x = y
            end if
         end do

         if (i /= n) then
             Root (i) = Root (n)
             Root (n) = x
         end if

      end do

      m = 1                                     ! elimination of duplicates
      do n = 2,nRoot
         x = Root (n)
         y = Root (m)
         z = abs (x - y)
         if (z > abs (x) * accuracy) then
             m = m + 1
             Root (m) = x
         end if
      end do

      nRoot = m

  end if

  if (printInfo) then
      do n = 1,nRoot
         write (*,'(a,es24.16)') ' Final quartic root    = ',Root (n)
      end do
      write (*,'(a)') ' ------------------------------------------------'
  end if
!
!
!     ...Transmit all roots.
!
!
  if (nRoot == 1) then

      Root1 = Root (1)

  else if (nRoot == 2) then

      Root1 = Root (1)
      Root2 = Root (2)

  else if (nRoot == 3) then

      Root1 = Root (1)
      Root2 = Root (2)
      Root3 = Root (3)

  else 

      Root1 = Root (1)
      Root2 = Root (2)
      Root3 = Root (3)
      Root4 = Root (4)

  end if
!
!
!     ...Ready!
!
!
  return
end subroutine ut_quarticRealRoots
