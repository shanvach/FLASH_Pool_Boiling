!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Integrate/save/op_BiggsPlanckIntegrate
!!
!! NAME
!!
!!  op_BiggsPlanckIntegrate
!!
!! SYNOPSIS
!!
!!  call op_BiggsPlanckIntegrate (real (in)  :: A1,
!!                                real (in)  :: A2,
!!                                real (in)  :: A3,
!!                                real (in)  :: A4,
!!                                real (in)  :: R,
!!                                real (in)  :: S,
!!                                real (in)  :: eBaseLowestExponent,
!!                                real (in)  :: eBaseLargestExponent,
!!                                real (in)  :: eBaseRescaleExponent,
!!                                real (out) :: ebaseScalingExponent,
!!                                real (out) :: BiggsPlanckIntegral,
!!                                real (out) :: PlanckIntegral)
!!
!! DESCRIPTION
!!
!!  This routine evaluates the following two kinds of integrals using the Romberg
!!  integration technique:
!!
!!                                          S
!!                                         /
!!                                        |
!!             Biggs-Planck Integral =    |  (A1/x + A2/x^2 + A3/x^3 + A4/x^4) * x^3 / (e^x - 1) dx
!!                                        |
!!                                       /
!!                                      R
!!
!!                                          S
!!                                         /
!!                                        |
!!                   Planck Integral =    |   x^3 / (e^x - 1) dx
!!                                        |
!!                                       /
!!                                      R
!!
!!  Several points need to be observed when evaluating these integrals via the Romberg discrete
!!  integration scheme. The problems come from the exponential function e^x present in both
!!  integrals, which for large x leads to overflow on the machine. Since the range of x can
!!  be very substantial, special measures have to be taken to get reasonable answers.
!!
!!  As the lower integration limit R gets larger, the exponential function in both integrals
!!  gets closer to the point of computational overflow. At a certain R value it is then advantageous
!!  to introduce the rescaling factor e^R, such that we have representable exponential values for a
!!  wide range in x. Both integrals are then multiplied by e^R and the resulting factor e^R/(e^x - 1)
!!  in the integrand is approximated as e^R/e^x = e^(R-x). Hence rescaling can only be applied when
!!  it is safe to ignore the -1 in (e^x - 1), meaning that this truncation does not lead to accuracy
!!  errors. It is hence dependent on the number of mantissa digits. The value of this limiting R for
!!  when the rescaling has to kick in is simply passsed as an argument and must have been determined
!!  elsewere.
!!
!!  The upper integration limit S might also be way beyond representability of e^S. The code has
!!  thus to decide upon a different upper integration limit U < S, beyond which all integrands are
!!  considered to be equal to zero. U is determined from the largest possible exponent (passed in argument).
!!  If the condition U < S is found, S is replaced by U in the above integrals.
!!
!! ARGUMENTS
!!
!!  A1                   : 1st Biggs expansion coefficient for the x^(-1) function
!!  A2                   : 2nd Biggs expansion coefficient for the x^(-2) function
!!  A3                   : 3rd Biggs expansion coefficient for the x^(-3) function
!!  A4                   : 4th Biggs expansion coefficient for the x^(-4) function
!!  R                    : lower integration limit
!!  S                    : upper integration limit
!!  eBaseLowestExponent  : the lowest -x for which e^(-x) is still representable on the machine
!!  eBaseLargestExponent : the largest x for which e^(+x) is still representable on the machine
!!  eBaseRescaleExponent : the limiting x in e^(+x) above which rescaling of integrals is done
!!  eBaseScalingExponent : the value R in e^R for the rescaling factor (if needed)
!!  BiggsPlanckIntegral  : the value of the Biggs-Planck integral
!!  PlanckIntegral       : the value of the Planck integral
!!
!!***

subroutine op_BiggsPlanckIntegrate (A1,A2,A3,A4,                                 &
                                    R,S,                                         &
                                    eBaseLowestExponent,                         &
                                    eBaseLargestExponent,                        &
                                    eBaseRescaleExponent,                        &
                                                           eBaseScalingExponent, &
                                                           BiggsPlanckIntegral,  &
                                                           PlanckIntegral)

  use Driver_interface,      ONLY : Driver_abortFlash
  use Opacity_dataIntegrate, ONLY : op_minRombergSteps,   &
                                    op_maxRombergSteps,   &
                                    op_RombergAccuracy,   &
                                    op_RombergIntegrals1, &
                                    op_RombergIntegrals2, &
                                    op_RombergRows1,      &
                                    op_RombergRows2

  implicit none

# include "Opacity.h"

  real, intent (in)  :: A1,A2,A3,A4
  real, intent (out) :: BiggsPlanckIntegral
  real, intent (in)  :: eBaseLowestExponent
  real, intent (in)  :: eBaseLargestExponent
  real, intent (in)  :: eBaseRescaleExponent
  real, intent (out) :: eBaseScalingExponent
  real, intent (out) :: PlanckIntegral
  real, intent (in)  :: R,S

  logical :: accurate
  logical :: convergedBiggsPlanck
  logical :: convergedPlanck
  logical :: rescaledIntegral

  integer :: k,m,n
  integer :: nIntEvalSteps
  integer :: nRombergStepsPlanck
  integer :: nRombergStepsBiggsPlanck

  real    :: accuracy
  real    :: exphh,expRh
  real    :: h,hh,Rh
  real    :: Intold,Intnew
  real    :: invExpR,invExpU,invExpUR
  real    :: invExph,invExphh,invExpRh
  real    :: p4m,invp4m
  real    :: range
  real    :: sumValueP,sumValueBP
  real    :: U

  real, parameter :: zero = 0.0
  real, parameter :: half = 0.5
  real, parameter :: one  = 1.0
  real, parameter :: four = 4.0
!
!
!   ...Set default convergence.
!
!
  convergedPlanck      = .false.
  convergedBiggsPlanck = .false.

  nRombergStepsPlanck      = op_maxRombergSteps
  nRombergStepsBiggsPlanck = op_maxRombergSteps
!
!
!   ...Set the rescaling factor and the upper integration limit U.
!
!
  if (R < eBaseRescaleExponent) then
      rescaledIntegral     = .false.
      eBaseScalingExponent = zero
  else
      rescaledIntegral     = .true.
      eBaseScalingExponent = R
  end if

  U = min (eBaseScalingExponent - eBaseLowestExponent , S )

  if (U <= R) then
      call Driver_abortFlash ('[op_BiggsPlanckIntegrate] ERROR: Integration limit error: U =< R')
  end if
!
!
!   ...Start Romberg integration on both integrals.
!
!
  write (*,*) ' R = ',R
  write (*,*) ' U = ',U
  write (*,*) ' rescaled = ',rescaledIntegral

  range = U - R

  if (rescaledIntegral) then

      invExpUR   = exp (R - U)
      sumValueP  =                R*R*R  + invExpUR *  U*U*U
      sumValueBP = (R*(A2+A1*R)+A3+A4/R) + invExpUR * (U*(A2+A1*U)+A3+A4/U)
  else
      invExpR    = one / (exp (R) - one)
      invExpU    = one / (exp (U) - one)
      sumValueP  = invExpR *  R*R*R                + invExpU *  U*U*U
      sumValueBP = invExpR * (R*(A2+A1*R)+A3+A4/R) + invExpU * (U*(A2+A1*U)+A3+A4/U)
  end if

  op_RombergRows1 (0,LOW) = half * range * sumValueP
  op_RombergRows2 (0,LOW) = half * range * sumValueBP

  op_RombergIntegrals1 (0) = op_RombergRows1 (0,LOW)
  op_RombergIntegrals2 (0) = op_RombergRows2 (0,LOW)
!
!
!   ...Romberg iterates on both integrals.
!
!
  h = range

  do n = 1,op_maxRombergSteps

     h  = h * half
     hh = h + h
     Rh = R + h

     nIntEvalSteps = 2 ** (n-1) - 1
!
!
!   ...Rescaled integral evaluation. Evaluate only those that are needed.
!
!
     if (rescaledIntegral) then

         invExph  = exp (-h)
         invExphh = invExph * invExph

         if (.not.convergedPlanck .and. .not.convergedBiggsPlanck) then

              sumValueP  = invExph *  Rh*Rh*Rh
              sumValueBP = invExph * (Rh*(A2+A1*Rh)+A3+A4/Rh)

              do k = 1,nIntEvalSteps
                 Rh         = Rh + hh
                 invExph    = invExph * invExphh
                 sumValueP  = sumValueP  + invExph *  Rh*Rh*Rh
                 sumValueBP = sumValueBP + invExph * (Rh*(A2+A1*Rh)+A3+A4/Rh)
              end do

         else if (.not.convergedPlanck) then

              sumValueP  = invExph *  Rh*Rh*Rh

              do k = 1,nIntEvalSteps
                 Rh         = Rh + hh
                 invExph    = invExph * invExphh
                 sumValueP  = sumValueP  + invExph *  Rh*Rh*Rh
              end do

         else if (.not.convergedBiggsPlanck) then

              sumValueBP = invExph * (Rh*(A2+A1*Rh)+A3+A4/Rh)

              do k = 1,nIntEvalSteps
                 Rh         = Rh + hh
                 invExph    = invExph * invExphh
                 sumValueBP = sumValueBP + invExph * (Rh*(A2+A1*Rh)+A3+A4/Rh)
              end do

         end if
!
!
!   ...Normal integral evaluation. Evaluate only those that are needed.
!
!
     else

         exphh = exp (hh)
         expRh = exp (Rh)

         if (.not.convergedPlanck .and. .not.convergedBiggsPlanck) then

              invExpRh   = one / (expRh - one)
              sumValueP  = invExpRh *  Rh*Rh*Rh
              sumValueBP = invExpRh * (Rh*(A2+A1*Rh)+A3+A4/Rh)

              do k = 1,nIntEvalSteps
                 Rh         = Rh + hh
                 expRh      = expRh * exphh
                 invExpRh   = one / (expRh - one)
                 sumValueP  = sumValueP  + invExpRh *  Rh*Rh*Rh
                 sumValueBP = sumValueBP + invExpRh * (Rh*(A2+A1*Rh)+A3+A4/Rh)
              end do

         else if (.not.convergedPlanck) then

              sumValueP = Rh*Rh*Rh / (expRh - one)

              do k = 1,nIntEvalSteps
                 Rh        = Rh + hh
                 expRh     = expRh * exphh
                 sumValueP = sumValueP  +  Rh*Rh*Rh / (expRh - one)
              end do

         else if (.not.convergedBiggsPlanck) then

              sumValueBP = (Rh*(A2+A1*Rh)+A3+A4/Rh) / (expRh - one)

              do k = 1,nIntEvalSteps
                 Rh         = Rh + hh
                 expRh      = expRh * exphh
                 sumValueBP = sumValueBP + (Rh*(A2+A1*Rh)+A3+A4/Rh) / (expRh - one)
              end do

         end if

     end if
!
!
!   ...Construction of new n-th Romberg row.
!
!
     if (.not.convergedBiggsPlanck) then

          op_RombergRows2 (0,HIGH) = half * op_RombergRows2 (0,LOW) + h * sumValueBP

          p4m = one
          do m = 1,n
             p4m    = p4m * four
             invp4m = one / (p4m - one)
             op_RombergRows2 (m,HIGH) = (p4m * op_RombergRows2 (m-1,HIGH) - op_RombergRows2 (m-1,LOW)) * invp4m
          end do

          op_RombergIntegrals2 (n) = op_RombergRows2 (n,HIGH)

     end if

     if (.not.convergedPlanck) then

          op_RombergRows1 (0,HIGH) = half * op_RombergRows1 (0,LOW) + h * sumValueP

          p4m = one
          do m = 1,n
             p4m    = p4m * four
             invp4m = one / (p4m - one)
             op_RombergRows1 (m,HIGH) = (p4m * op_RombergRows1 (m-1,HIGH) - op_RombergRows1 (m-1,LOW)) * invp4m
          end do

          op_RombergIntegrals1 (n) = op_RombergRows1 (n,HIGH)

     end if
!
!
!   ...check accuracy of the integrals at the n-th Romberg step. The set of the latest
!      four integrals are analyzed in terms of relative accuracy obtained at each of the
!      last four steps. If all three relative accuracies are below the specified accuracy
!      limit, the integral(s) is(are) considered to be converged.
!
!
     if (n > op_minRombergSteps) then

         if (.not.convergedPlanck) then
             accurate  = .true.
             do m = n-2,n
                Intold    = op_RombergIntegrals1 (m-1)
                Intnew    = op_RombergIntegrals1 (m)
                accuracy  = abs ((Intnew - Intold) / Intnew)
                accurate  = accurate  .and. (accuracy  < op_RombergAccuracy)
             end do
             convergedPlanck = accurate
         end if

         if (.not.convergedBiggsPlanck) then
             accurate  = .true.
             do m = n-2,n
                Intold    = op_RombergIntegrals2 (m-1)
                Intnew    = op_RombergIntegrals2 (m)
                accuracy  = abs ((Intnew - Intold) / Intnew)
                accurate  = accurate  .and. (accuracy  < op_RombergAccuracy)
             end do
             convergedBiggsPlanck = accurate
         end if

     end if

     if (convergedPlanck .and. convergedBiggsPlanck) then
         nRombergStepsPlanck      = min (n,nRombergStepsPlanck)
         nRombergStepsBiggsPlanck = min (n,nRombergStepsBiggsPlanck)
         exit
     else if (convergedPlanck) then
         nRombergStepsPlanck       = min (n,nRombergStepsPlanck)
         op_RombergRows2 (0:n,LOW) = op_RombergRows2 (0:n,HIGH)
     else if (convergedBiggsPlanck) then
         nRombergStepsBiggsPlanck  = min (n,nRombergStepsBiggsPlanck)
         op_RombergRows1 (0:n,LOW) = op_RombergRows1 (0:n,HIGH)
     else
         op_RombergRows1 (0:n,LOW) = op_RombergRows1 (0:n,HIGH)
         op_RombergRows2 (0:n,LOW) = op_RombergRows2 (0:n,HIGH)
     end if        
    
  end do

  if (.not.convergedPlanck) then
       Intold    = op_RombergIntegrals1 (op_maxRombergSteps-1)
       Intnew    = op_RombergIntegrals1 (op_maxRombergSteps)
       accuracy  = abs ((Intnew - Intold) / Intnew)
       write (*,*) ' Planck integral not converged after ',op_maxRombergSteps,' Romberg steps! '
       write (*,*) ' Highest Planck integral relative accuracy = ',accuracy
       write (*,*) '                           Planck Integral = ',Intnew
  else
       PlanckIntegral = op_RombergIntegrals1 (nRombergStepsPlanck)
       write (*,*) ' # of       Planck Romberg steps = ',nRombergStepsPlanck
       write (*,*) '                 Planck Integral = ',PlanckIntegral
  end if

  if (.not.convergedBiggsPlanck) then
       Intold    = op_RombergIntegrals2 (op_maxRombergSteps-1)
       Intnew    = op_RombergIntegrals2 (op_maxRombergSteps)
       accuracy  = abs ((Intnew - Intold) / Intnew)
       write (*,*) ' Biggs-Planck integral not converged after ',op_maxRombergSteps,' Romberg steps! '
       write (*,*) ' Highest Biggs-Planck integral relative accuracy = ',accuracy
       write (*,*) '                           Biggs-Planck Integral = ',Intnew
  else
       BiggsPlanckIntegral = op_RombergIntegrals2 (nRombergStepsBiggsPlanck)
       write (*,*) ' # of Biggs Planck Romberg steps = ',nRombergStepsBiggsPlanck
       write (*,*) '           Biggs-Planck Integral = ',BiggsPlanckIntegral
  end if
!
!
!   ...Ready! 
!
!
  return
end subroutine op_BiggsPlanckIntegrate
