!!****if* source/flashUtilities/1dRootfind/ut_bracketed_1d_rootfind
!!
!! NAME
!!
!!  bracketed_1d_rootfind
!!
!! SYNOPSIS
!!
!!  bracketed_1d_rootfind(fun, tgt_val, x_l, x_r, f_l, f_r,
!!                           tolerance, x_soln, error, deriv, ncalls, niter,
!!                           absolute_convergence, root_convergence)
!!
!! DESCRIPTION
!!
!! This routine locates the root x_soln of the equation
!! fun(x_soln)=tgt_val, to relative precision rel_precision. The
!! user supplied bounds x_l and x_r must bracket the root, that is, they
!! must satisfy (fun(x_l) - tgt_val)*(fun(x_r) - tgt_val) < 0.  Failure
!! to satisfy this input condition results in immediate return with
!! error=1.  The condition is actually tested using f_l and f_r, assuming
!! the calling routine set these to f_l=fun(x_l), f_r=fun(x_r).
!!
!! If the optional argument deriv (a function of one real argument,
!! like fun) is supplied, then the algorithm is a hybrid Newton-Raphson/
!! bisection minimization, clamped to stay within the (shrinking)
!! bracketed interval.  If deriv is not present, the algorithm is
!! straight bisection.
!!
!! Note that while possible, it is not a great idea to use a finite-
!! difference derivative algorithm for deriv.  As discussed in
!! Numerical Recipes, doing so reduces the effective order of
!! convergence (per function evaluation).  Worse, the differencing
!! interval must be greater than roundoff by the relative precision,
!! but if it is too large then the convergence is no better than
!! linear anyway, so you might as well just use bisection.  Generally,
!! Newton-Raphson is for problems where the derivative is computable
!! directly.
!!
!! The root is returned in x_soln, and the total number of function
!! calls (including to deriv, if any) is returned in the optional
!! argument ncalls.  The maximum number of bisection iterations is 40.
!! Root searches exceeding this maximum result in a return with error=2.
!! Convergence criteria are controlled by the arguments 'tolerance',
!! 'absolute_convergence', and 'root_convergence', as described below.
!!
!! On success, error=0.
!!
!! ARGUMENTS:
!!
!! fun : Function to be rooted
!!
!! tgt_val : Value of fun at the desired root
!!
!! x_l : Left (i.e. lower) bound of real interval bracketing the root
!!
!! x_r : Right (i.e. higher) bound of said interval
!!
!! f_l : Value of fun(x_l)
!!
!! f_r : Value of fun(x_r)
!!
!! tolerance :     Precision of root location.  By default, convergence
!!                 is declared if abs(fun(x_soln) - tgt_val)/tgt_val < tolerance.
!!                 This behavior is different if the optional argument
!!                 absolute_convergence is present.
!!
!! x_soln : The root.  If the routine fails, this will be set to
!!          huge(x_soln), to get the attention of someone in authority
!!          even if nobody is checking the value of error (so there).
!!
!! error : On success, error=0. If x_l and x_r do not bracket the root
!!         on input, error=1.  If the maximum number of iterations is
!!         exceeded, error=2.
!!
!! deriv : Optional, function returning the derivative of fun
!!
!! ncalls: Number of calls to fun and to derivs (optional)
!!
!! niter: Number of bisection iterations (optional)
!!
!! absolute_convergence: Optional, if present the convergence criterion
!!                       becomes abs(fun(x_soln) - tgt_val) < tolerance
!!                       (if .not.present(root_convergence) ) or
!!                       abs(x_soln - previous_x_soln) < tolerance
!!                       (otherwise).  By default, the convergence
!!                       criterion is the relative version of these.
!!
!! root_convergence: Optional, if present the convergence criterion is
!!                   according to the root value, rather than to the
!!                   function value.
!!
!!***

subroutine bracketed_1d_rootfind(fun, tgt_val, x_l, x_r, f_l, f_r, &
                                 tolerance, x_soln, error, deriv, ncalls, &
                                 niter, absolute_convergence, root_convergence)

implicit none

interface

  real function fun(x)
    real, intent(in) :: x
  end function fun

  real function deriv(x)
    real, intent(in) :: x
  end function deriv

end interface

optional :: deriv

integer, intent(in), optional :: absolute_convergence, root_convergence
real, intent(in) :: tgt_val, tolerance
real, intent(inout) :: x_l, x_r, f_l, f_r
real, intent(out) :: x_soln
integer, intent(out) :: error
integer, optional, intent(out) :: ncalls, niter

integer, parameter :: maxiter_b = 40, maxiter_nr = 20
real, parameter :: maxdx_nr = 0.1
integer :: iter_b, iter_nr, nc
logical :: do_nr, nr_step, converged, relative, by_funval
real :: x_curr, x_b, x_nr, delx_b, delx_nr, f_curr, df_curr, x_prev
real :: v1, v2

  if ( (f_l - tgt_val) * (f_r - tgt_val) .gt. 0.0 ) then
    error = 1
    x_soln = huge(x_soln)
    if ( present(ncalls) ) ncalls = 0
    return
  endif
  do_nr = present(deriv)
  relative = .not. present(absolute_convergence)
  by_funval = .not. present(root_convergence)
  
  x_curr = 0.5 * (x_l + x_r)
  f_curr = fun(x_curr)
  nc = 1
  converged = .false.
  iter_nr = 0
  iter_b = 0

  do
  
    x_prev = x_curr
    if ( (f_curr - tgt_val) * (f_l - tgt_val) .lt. 0.0 ) then
      x_b = 0.5 * (x_l + x_curr)
      x_r = x_curr
    else
      x_b = 0.5 * (x_r + x_curr)
      x_l = x_curr
    endif
    delx_b = x_b - x_curr
    
    nr_step = do_nr
    if ( nr_step ) then

      df_curr = deriv(x_curr)
      nc = nc + 1
      nr_step = nr_step .and. ( abs(df_curr) .gt. tiny(df_curr) * abs(f_curr - tgt_val) ) 
                                                        !! If the derivative is tiny, don't bother with NR
    endif
    
    if ( nr_step ) then
      delx_nr = ( tgt_val - f_curr ) / df_curr
      x_nr = x_curr + delx_nr
      !! Newton-Raphson step ONLY if:
      nr_step = nr_step .and. ( (delx_b*delx_nr) .gt. 0.0 ) & !! Goes in same direction as bisection step;
                        .and. ( abs(delx_nr) .lt. maxdx_nr * abs(x_curr) ) & !! Relative step size does not exceed limit;
                        .and. ( x_nr .gt. x_l ) .and. (x_nr .lt. x_r ) & !! Step doesn't leave bracketed interval;
                        .and. ( iter_nr .le. maxiter_nr ) !! There haven't been too many NR steps already.
    endif
    
    if ( nr_step ) then
      x_curr = x_nr
      iter_nr = iter_nr + 1
    else
      x_curr = x_b
      iter_b = iter_b + 1
      iter_nr = 0
    endif
    f_curr = fun(x_curr)
    nc = nc + 1
    
    if ( by_funval ) then
      v1 = abs(f_curr - tgt_val)
      v2 = abs(tgt_val)
    else
      v1 = abs(x_curr - x_prev)
      v2 = 0.5*abs(x_curr + x_prev)
    endif

    if (relative) then
      converged = ( v1 .lt. tolerance * v2 )
    else
      converged = ( v1 .lt. tolerance )
    endif
    
    if ( converged .or. iter_b .gt. maxiter_b ) exit
  
  end do
  
  if ( converged ) then
    error = 0
    x_soln = x_curr
  else
    error = 2
    x_soln = huge(x_soln)
  endif
  if ( present(ncalls) ) ncalls = nc
  if ( present(niter) ) niter = iter_b
  
  return
  
end subroutine bracketed_1d_rootfind
