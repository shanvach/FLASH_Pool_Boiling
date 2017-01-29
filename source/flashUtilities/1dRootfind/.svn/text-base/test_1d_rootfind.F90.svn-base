  program test_1d_rootfind
  
  use ut_bracketed_1d_rootfindinterface

  implicit none

  real :: fn_cont, fn_disc
  external  :: fn_cont, fn_disc

  real :: tol = 1.0E-07, tgt_val = 1.0, x_l = 1.0E-02, x_r = 0.99
  real :: xll, xrr, fl, fr, x_soln, f_soln
  integer :: ni, nc, error

  write(*,*) 'Equation to be rooted: tan(pi*x/2) - pi*x/2 = ', tgt_val
  write(*,*) 'Tolerance: ', tol
  write(*,*) 'Initial root bracket: x_l = ', x_l, ', x_r = ', x_r
  write(*,*) 

  xll = x_l
  xrr = x_r
  fl = fn_cont(xll)
  fr = fn_cont(xrr)
  call cont_rel_byfun_wder(tgt_val, tol, xll, xrr, fl, fr, &
                           x_soln, f_soln, error, ni, nc)
  write(*,*) 'Continuous Function, Relative Convergence by Function Value with Derivatives:'
  if ( error .ne. 0 ) then
    write(*,*) 'bracketed_1d_rootfind returned error code ', error
    stop
  endif
  f_soln = fn_cont(x_soln)
  write(*,'(a, 1pE15.7, a, 1pE15.7)') 'Root = ', x_soln, ', Function = ', f_soln
  write(*,*) '# of Bisection Calls = ', ni, ' # of Function Calls = ', nc
  write(*,*)

  xll = x_l
  xrr = x_r
  fl = fn_cont(xll)
  fr = fn_cont(xrr)
  call cont_rel_byfun_woutder(tgt_val, tol, xll, xrr, fl, fr, &
                              x_soln, f_soln, error, ni, nc)
  write(*,*) 'Continuous Function, Relative Convergence by Function Value without Derivatives:'
  if ( error .ne. 0 ) then
    write(*,*) 'bracketed_1d_rootfind returned error code ', error
    stop
  endif
  f_soln = fn_cont(x_soln)
  write(*,'(a, 1pE15.7, a, 1pE15.7)') 'Root = ', x_soln, ', Function = ', f_soln
  write(*,*) '# of Bisection Calls = ', ni, ' # of Function Calls = ', nc
  write(*,*) 

  xll = x_l
  xrr = x_r
  fl = fn_cont(xll)
  fr = fn_cont(xrr)
  call cont_rel_byroot_wder(tgt_val, tol, xll, xrr, fl, fr, &
                           x_soln, f_soln, error, ni, nc)
  write(*,*) 'Continuous Function, Relative Convergence by Root Value with Derivatives:'
  if ( error .ne. 0 ) then
    write(*,*) 'bracketed_1d_rootfind returned error code ', error
    stop
  endif
  f_soln = fn_cont(x_soln)
  write(*,'(a, 1pE15.7, a, 1pE15.7)') 'Root = ', x_soln, ', Function = ', f_soln
  write(*,*) '# of Bisection Calls = ', ni, ' # of Function Calls = ', nc
  write(*,*)

  xll = x_l
  xrr = x_r
  fl = fn_cont(xll)
  fr = fn_cont(xrr)
  call cont_rel_byroot_woutder(tgt_val, tol, xll, xrr, fl, fr, &
                           x_soln, f_soln, error, ni, nc)
  write(*,*) 'Continuous Function, Relative Convergence by Root Value without Derivatives:'
  if ( error .ne. 0 ) then
    write(*,*) 'bracketed_1d_rootfind returned error code ', error
    stop
  endif
  f_soln = fn_cont(x_soln)
  write(*,'(a, 1pE15.7, a, 1pE15.7)') 'Root = ', x_soln, ', Function = ', f_soln
  write(*,*) '# of Bisection Calls = ', ni, ' # of Function Calls = ', nc
  write(*,*)

  xll = x_l
  xrr = x_r
  fl = fn_cont(xll)
  fr = fn_cont(xrr)
  call cont_abs_byfun_wder(tgt_val, tol, xll, xrr, fl, fr, &
                           x_soln, f_soln, error, ni, nc)
  write(*,*) 'Continuous Function, Absolute Convergence by Function Value with Derivatives:'
  if ( error .ne. 0 ) then
    write(*,*) 'bracketed_1d_rootfind returned error code ', error
    stop
  endif
  f_soln = fn_cont(x_soln)
  write(*,'(a, 1pE15.7, a, 1pE15.7)') 'Root = ', x_soln, ', Function = ', f_soln
  write(*,*) '# of Bisection Calls = ', ni, ' # of Function Calls = ', nc
  write(*,*)

  xll = x_l
  xrr = x_r
  fl = fn_cont(xll)
  fr = fn_cont(xrr)
  call cont_abs_byfun_woutder(tgt_val, tol, xll, xrr, fl, fr, &
                           x_soln, f_soln, error, ni, nc)
  write(*,*) 'Continuous Function, Absolute Convergence by Function Value without Derivatives:'
  if ( error .ne. 0 ) then
    write(*,*) 'bracketed_1d_rootfind returned error code ', error
    stop
  endif
  f_soln = fn_cont(x_soln)
  write(*,'(a, 1pE15.7, a, 1pE15.7)') 'Root = ', x_soln, ', Function = ', f_soln
  write(*,*) '# of Bisection Calls = ', ni, ' # of Function Calls = ', nc
  write(*,*)

  xll = x_l
  xrr = x_r
  fl = fn_cont(xll)
  fr = fn_cont(xrr)
  call cont_abs_byroot_wder(tgt_val, tol, xll, xrr, fl, fr, &
                           x_soln, f_soln, error, ni, nc)
  write(*,*) 'Continuous Function, Absolute Convergence by Root Value with Derivatives:'
  if ( error .ne. 0 ) then
    write(*,*) 'bracketed_1d_rootfind returned error code ', error
    stop
  endif
  f_soln = fn_cont(x_soln)
  write(*,'(a, 1pE15.7, a, 1pE15.7)') 'Root = ', x_soln, ', Function = ', f_soln
  write(*,*) '# of Bisection Calls = ', ni, ' # of Function Calls = ', nc
  write(*,*)

  xll = x_l
  xrr = x_r
  fl = fn_cont(xll)
  fr = fn_cont(xrr)
  call cont_abs_byroot_woutder(tgt_val, tol, xll, xrr, fl, fr, &
                           x_soln, f_soln, error, ni, nc)
  write(*,*) 'Continuous Function, Absolute Convergence by Root Value without Derivatives:'
  if ( error .ne. 0 ) then
    write(*,*) 'bracketed_1d_rootfind returned error code ', error
    stop
  endif
  f_soln = fn_cont(x_soln)
  write(*,'(a, 1pE15.7, a, 1pE15.7)') 'Root = ', x_soln, ', Function = ', f_soln
  write(*,*) '# of Bisection Calls = ', ni, ' # of Function Calls = ', nc
  write(*,*)

!!!      !!!     !!!       !!!

  xll = x_l
  xrr = x_r
  fl = fn_disc(xll)
  fr = fn_disc(xrr)
  call disc_rel_byfun_wder(tgt_val, tol, xll, xrr, fl, fr, &
                           x_soln, f_soln, error, ni, nc)
  write(*,*) 'Discretized Function, Relative Convergence by Function Value with Derivatives:'
  if ( error .ne. 0 ) then
    write(*,*) 'bracketed_1d_rootfind returned error code ', error
    stop
  endif
  f_soln = fn_disc(x_soln)
  write(*,'(a, 1pE15.7, a, 1pE15.7)') 'Root = ', x_soln, ', Function = ', f_soln
  write(*,*) '# of Bisection Calls = ', ni, ' # of Function Calls = ', nc
  write(*,*)

  xll = x_l
  xrr = x_r
  fl = fn_disc(xll)
  fr = fn_disc(xrr)
  call disc_rel_byfun_woutder(tgt_val, tol, xll, xrr, fl, fr, &
                              x_soln, f_soln, error, ni, nc)
  write(*,*) 'Discretized Function, Relative Convergence by Function Value without Derivatives:'
  if ( error .ne. 0 ) then
    write(*,*) 'bracketed_1d_rootfind returned error code ', error
    stop
  endif
  f_soln = fn_disc(x_soln)
  write(*,'(a, 1pE15.7, a, 1pE15.7)') 'Root = ', x_soln, ', Function = ', f_soln
  write(*,*) '# of Bisection Calls = ', ni, ' # of Function Calls = ', nc
  write(*,*) 

  xll = x_l
  xrr = x_r
  fl = fn_disc(xll)
  fr = fn_disc(xrr)
  call disc_rel_byroot_wder(tgt_val, tol, xll, xrr, fl, fr, &
                           x_soln, f_soln, error, ni, nc)
  write(*,*) 'Discretized Function, Relative Convergence by Root Value with Derivatives:'
  if ( error .ne. 0 ) then
    write(*,*) 'bracketed_1d_rootfind returned error code ', error
    stop
  endif
  f_soln = fn_disc(x_soln)
  write(*,'(a, 1pE15.7, a, 1pE15.7)') 'Root = ', x_soln, ', Function = ', f_soln
  write(*,*) '# of Bisection Calls = ', ni, ' # of Function Calls = ', nc
  write(*,*)

  xll = x_l
  xrr = x_r
  fl = fn_disc(xll)
  fr = fn_disc(xrr)
  call disc_rel_byroot_woutder(tgt_val, tol, xll, xrr, fl, fr, &
                           x_soln, f_soln, error, ni, nc)
  write(*,*) 'Discretized Function, Relative Convergence by Root Value without Derivatives:'
  if ( error .ne. 0 ) then
    write(*,*) 'bracketed_1d_rootfind returned error code ', error
    stop
  endif
  f_soln = fn_disc(x_soln)
  write(*,'(a, 1pE15.7, a, 1pE15.7)') 'Root = ', x_soln, ', Function = ', f_soln
  write(*,*) '# of Bisection Calls = ', ni, ' # of Function Calls = ', nc
  write(*,*)

  xll = x_l
  xrr = x_r
  fl = fn_disc(xll)
  fr = fn_disc(xrr)
  call disc_abs_byfun_wder(tgt_val, tol, xll, xrr, fl, fr, &
                           x_soln, f_soln, error, ni, nc)
  write(*,*) 'Discretized Function, Absolute Convergence by Function Value with Derivatives:'
  if ( error .ne. 0 ) then
    write(*,*) 'bracketed_1d_rootfind returned error code ', error
    stop
  endif
  f_soln = fn_disc(x_soln)
  write(*,'(a, 1pE15.7, a, 1pE15.7)') 'Root = ', x_soln, ', Function = ', f_soln
  write(*,*) '# of Bisection Calls = ', ni, ' # of Function Calls = ', nc
  write(*,*)

  xll = x_l
  xrr = x_r
  fl = fn_disc(xll)
  fr = fn_disc(xrr)
  call disc_abs_byfun_woutder(tgt_val, tol, xll, xrr, fl, fr, &
                           x_soln, f_soln, error, ni, nc)
  write(*,*) 'Discretized Function, Absolute Convergence by Function Value without Derivatives:'
  if ( error .ne. 0 ) then
    write(*,*) 'bracketed_1d_rootfind returned error code ', error
    stop
  endif
  f_soln = fn_disc(x_soln)
  write(*,'(a, 1pE15.7, a, 1pE15.7)') 'Root = ', x_soln, ', Function = ', f_soln
  write(*,*) '# of Bisection Calls = ', ni, ' # of Function Calls = ', nc
  write(*,*)

  xll = x_l
  xrr = x_r
  fl = fn_disc(xll)
  fr = fn_disc(xrr)
  call disc_abs_byroot_wder(tgt_val, tol, xll, xrr, fl, fr, &
                           x_soln, f_soln, error, ni, nc)
  write(*,*) 'Discretized Function, Absolute Convergence by Root Value with Derivatives:'
  if ( error .ne. 0 ) then
    write(*,*) 'bracketed_1d_rootfind returned error code ', error
    stop
  endif
  f_soln = fn_disc(x_soln)
  write(*,'(a, 1pE15.7, a, 1pE15.7)') 'Root = ', x_soln, ', Function = ', f_soln
  write(*,*) '# of Bisection Calls = ', ni, ' # of Function Calls = ', nc
  write(*,*)

  xll = x_l
  xrr = x_r
  fl = fn_disc(xll)
  fr = fn_disc(xrr)
  call disc_abs_byroot_woutder(tgt_val, tol, xll, xrr, fl, fr, &
                           x_soln, f_soln, error, ni, nc)
  write(*,*) 'Discretized Function, Absolute Convergence by Root Value without Derivatives:'
  if ( error .ne. 0 ) then
    write(*,*) 'bracketed_1d_rootfind returned error code ', error
    stop
  endif
  f_soln = fn_disc(x_soln)
  write(*,'(a, 1pE15.7, a, 1pE15.7)') 'Root = ', x_soln, ', Function = ', f_soln
  write(*,*) '# of Bisection Calls = ', ni, ' # of Function Calls = ', nc
  write(*,*)

  end program test_1d_rootfind
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cont_rel_byfun_wder(tgt_val, tol, xl, xr, fl, fr, &
                                 x_soln, f_soln, error, ni, nc)

    use ut_bracketed_1d_rootfindinterface

    implicit none
    real, intent(in) :: tgt_val, tol
  real, intent(inout) :: xl, xr, fl, fr
    real, intent(out) :: x_soln, f_soln
    integer, intent(out) :: error, nc, ni
    real :: fn_cont, dfn_cont
    external  :: fn_cont, dfn_cont
    
    call bracketed_1d_rootfind(fn_cont, tgt_val, xl, xr, fl, fr, &
                               tol, x_soln, error, &
                               deriv=dfn_cont, &
                               ncalls=nc, &
                               niter=ni)
 
  end subroutine cont_rel_byfun_wder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cont_rel_byfun_woutder(tgt_val, tol, xl, xr, fl, fr, &
                                 x_soln, f_soln, error, ni, nc)

    use ut_bracketed_1d_rootfindinterface

    implicit none
    real, intent(in) :: tgt_val, tol
    real, intent(inout) :: xl, xr, fl, fr
    real, intent(out) :: x_soln, f_soln
    integer, intent(out) :: error, nc, ni
    real :: fn_cont
    external  :: fn_cont
 
    call bracketed_1d_rootfind(fn_cont, tgt_val, xl, xr, fl, fr, &
                               tol, x_soln, error,  &
                               ncalls=nc, &
                               niter=ni)
 
  end subroutine cont_rel_byfun_woutder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cont_rel_byroot_wder(tgt_val, tol, xl, xr, fl, fr, &
                                 x_soln, f_soln, error, ni, nc)

    use ut_bracketed_1d_rootfindinterface

    implicit none
    real, intent(in) :: tgt_val, tol
    real, intent(inout) :: xl, xr, fl, fr
    real, intent(out) :: x_soln, f_soln
    integer, intent(out) :: error, nc, ni
    real :: fn_cont, dfn_cont
    external  :: fn_cont, dfn_cont
    
    call bracketed_1d_rootfind(fn_cont, tgt_val, xl, xr, fl, fr, &
                               tol, x_soln, error,  &
                               deriv=dfn_cont, &
                               ncalls=nc, &
                               niter=ni, &
                               root_convergence=1)
 
  end subroutine cont_rel_byroot_wder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cont_rel_byroot_woutder(tgt_val, tol, xl, xr, fl, fr, &
                                 x_soln, f_soln, error, ni, nc)

    use ut_bracketed_1d_rootfindinterface

    implicit none
    real, intent(in) :: tgt_val, tol
    real, intent(inout) :: xl, xr, fl, fr
    real, intent(out) :: x_soln, f_soln
    integer, intent(out) :: error, nc, ni
    real :: fn_cont
    external  :: fn_cont
    
    call bracketed_1d_rootfind(fn_cont, tgt_val, xl, xr, fl, fr, &
                               tol, x_soln, error,  &
                               ncalls=nc, &
                               niter=ni, &
                               root_convergence=1)
 
  end subroutine cont_rel_byroot_woutder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cont_abs_byfun_wder(tgt_val, tol, xl, xr, fl, fr, &
                                 x_soln, f_soln, error, ni, nc)

    use ut_bracketed_1d_rootfindinterface

    implicit none
    real, intent(in) :: tgt_val, tol
    real, intent(inout) :: xl, xr, fl, fr
    real, intent(out) :: x_soln, f_soln
    integer, intent(out) :: error, nc, ni
    real :: fn_cont, dfn_cont
    external  :: fn_cont, dfn_cont
    
    call bracketed_1d_rootfind(fn_cont, tgt_val, xl, xr, fl, fr, &
                               tol, x_soln, error,  &
                               deriv=dfn_cont, &
                               ncalls=nc, &
                               niter=ni, &
                               absolute_convergence=1)
 
  end subroutine cont_abs_byfun_wder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cont_abs_byfun_woutder(tgt_val, tol, xl, xr, fl, fr, &
                                 x_soln, f_soln, error, ni, nc)

    use ut_bracketed_1d_rootfindinterface

    implicit none
    real, intent(in) :: tgt_val, tol
    real, intent(inout) :: xl, xr, fl, fr
    real, intent(out) :: x_soln, f_soln
    integer, intent(out) :: error, nc, ni
    real :: fn_cont
    external  :: fn_cont
    
    call bracketed_1d_rootfind(fn_cont, tgt_val, xl, xr, fl, fr, &
                               tol, x_soln, error,  &
                               ncalls=nc, &
                               niter=ni, &
                               absolute_convergence=1)
 
  end subroutine cont_abs_byfun_woutder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cont_abs_byroot_wder(tgt_val, tol, xl, xr, fl, fr, &
                                 x_soln, f_soln, error, ni, nc)

    use ut_bracketed_1d_rootfindinterface

    implicit none
    real, intent(in) :: tgt_val, tol
    real, intent(inout) :: xl, xr, fl, fr
    real, intent(out) :: x_soln, f_soln
    integer, intent(out) :: error, nc, ni
    real :: fn_cont, dfn_cont
    external  :: fn_cont, dfn_cont
    
    call bracketed_1d_rootfind(fn_cont, tgt_val, xl, xr, fl, fr, &
                               tol, x_soln, error,  &
                               deriv=dfn_cont, &
                               ncalls=nc, &
                               niter=ni, &
                               absolute_convergence=1, &
                               root_convergence=1)
 
  end subroutine cont_abs_byroot_wder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cont_abs_byroot_woutder(tgt_val, tol, xl, xr, fl, fr, &
                                 x_soln, f_soln, error, ni, nc)

    use ut_bracketed_1d_rootfindinterface

    implicit none
    real, intent(in) :: tgt_val, tol
    real, intent(inout) :: xl, xr, fl, fr
    real, intent(out) :: x_soln, f_soln
    integer, intent(out) :: error, nc, ni
    real :: fn_cont
    external  :: fn_cont
    
    call bracketed_1d_rootfind(fn_cont, tgt_val, xl, xr, fl, fr, &
                               tol, x_soln, error,  &
                               ncalls=nc, &
                               niter=ni, &
                               absolute_convergence=1, &
                               root_convergence=1)

   end subroutine cont_abs_byroot_woutder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine disc_rel_byfun_wder(tgt_val, tol, xl, xr, fl, fr, &
                                 x_soln, f_soln, error, ni, nc)

    use ut_bracketed_1d_rootfindinterface

    implicit none
    real, intent(in) :: tgt_val, tol
    real, intent(inout) :: xl, xr, fl, fr
    real, intent(out) :: x_soln, f_soln
    integer, intent(out) :: error, nc, ni
    real :: fn_disc, dfn_disc
    external  :: fn_disc, dfn_disc
    
    call bracketed_1d_rootfind(fn_disc, tgt_val, xl, xr, fl, fr, &
                               tol, x_soln, error,  &
                               deriv=dfn_disc, &
                               ncalls=nc, &
                               niter=ni)
 
  end subroutine disc_rel_byfun_wder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine disc_rel_byfun_woutder(tgt_val, tol, xl, xr, fl, fr, &
                                 x_soln, f_soln, error, ni, nc)

    use ut_bracketed_1d_rootfindinterface

    implicit none
    real, intent(in) :: tgt_val, tol
    real, intent(inout) :: xl, xr, fl, fr
    real, intent(out) :: x_soln, f_soln
    integer, intent(out) :: error, nc, ni
    real :: fn_disc
    external  :: fn_disc
 
    call bracketed_1d_rootfind(fn_disc, tgt_val, xl, xr, fl, fr, &
                               tol, x_soln, error,  &
                               ncalls=nc, &
                               niter=ni)
 
  end subroutine disc_rel_byfun_woutder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine disc_rel_byroot_wder(tgt_val, tol, xl, xr, fl, fr, &
                                 x_soln, f_soln, error, ni, nc)

    use ut_bracketed_1d_rootfindinterface

    implicit none
    real, intent(in) :: tgt_val, tol
    real, intent(inout) :: xl, xr, fl, fr
    real, intent(out) :: x_soln, f_soln
    integer, intent(out) :: error, nc, ni
    real :: fn_disc, dfn_disc
    external  :: fn_disc, dfn_disc
    
    call bracketed_1d_rootfind(fn_disc, tgt_val, xl, xr, fl, fr, &
                               tol, x_soln, error,  &
                               deriv= fn_disc, &
                               ncalls=nc, &
                               niter=ni, &
                               root_convergence=1)
 
  end subroutine disc_rel_byroot_wder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine disc_rel_byroot_woutder(tgt_val, tol, xl, xr, fl, fr, &
                                 x_soln, f_soln, error, ni, nc)

    use ut_bracketed_1d_rootfindinterface

    implicit none
    real, intent(in) :: tgt_val, tol
    real, intent(inout) :: xl, xr, fl, fr
    real, intent(out) :: x_soln, f_soln
    integer, intent(out) :: error, nc, ni
    real :: fn_disc
    external  :: fn_disc
    
    call bracketed_1d_rootfind(fn_disc, tgt_val, xl, xr, fl, fr, &
                               tol, x_soln, error,  &
                               ncalls=nc, &
                               niter=ni, &
                               root_convergence=1)
 
  end subroutine disc_rel_byroot_woutder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine disc_abs_byfun_wder(tgt_val, tol, xl, xr, fl, fr, &
                                 x_soln, f_soln, error, ni, nc)

    use ut_bracketed_1d_rootfindinterface

    implicit none
    real, intent(in) :: tgt_val, tol
    real, intent(inout) :: xl, xr, fl, fr
    real, intent(out) :: x_soln, f_soln
    integer, intent(out) :: error, nc, ni
    real :: fn_disc, dfn_disc
    external  :: fn_disc, dfn_disc
    
    call bracketed_1d_rootfind(fn_disc, tgt_val, xl, xr, fl, fr, &
                               tol, x_soln, error,  &
                               deriv=dfn_disc, &
                               ncalls=nc, &
                               niter=ni, &
                               absolute_convergence=1)
 
  end subroutine disc_abs_byfun_wder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine disc_abs_byfun_woutder(tgt_val, tol, xl, xr, fl, fr, &
                                 x_soln, f_soln, error, ni, nc)

    use ut_bracketed_1d_rootfindinterface

    implicit none
    real, intent(in) :: tgt_val, tol
    real, intent(inout) :: xl, xr, fl, fr
    real, intent(out) :: x_soln, f_soln
    integer, intent(out) :: error, nc, ni
    real :: fn_disc
    external  :: fn_disc
    
    call bracketed_1d_rootfind(fn_disc, tgt_val, xl, xr, fl, fr, &
                               tol, x_soln, error,  &
                               ncalls=nc, &
                               niter=ni, &
                               absolute_convergence=1)
 
  end subroutine disc_abs_byfun_woutder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine disc_abs_byroot_wder(tgt_val, tol, xl, xr, fl, fr, &
                                 x_soln, f_soln, error, ni, nc)

    use ut_bracketed_1d_rootfindinterface

    implicit none
    real, intent(in) :: tgt_val, tol
    real, intent(inout) :: xl, xr, fl, fr
    real, intent(out) :: x_soln, f_soln
    integer, intent(out) :: error, nc, ni
    real :: fn_disc, dfn_disc
    external  :: fn_disc, dfn_disc
    
    call bracketed_1d_rootfind(fn_disc, tgt_val, xl, xr, fl, fr, &
                               tol, x_soln, error,  &
                               deriv=dfn_disc, &
                               ncalls=nc, &
                               niter=ni, &
                               absolute_convergence=1, &
                               root_convergence=1)
 
  end subroutine disc_abs_byroot_wder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine disc_abs_byroot_woutder(tgt_val, tol, xl, xr, fl, fr, &
                                 x_soln, f_soln, error, ni, nc)

    use ut_bracketed_1d_rootfindinterface

    implicit none
    real, intent(in) :: tgt_val, tol
    real, intent(inout) :: xl, xr, fl, fr
    real, intent(out) :: x_soln, f_soln
    integer, intent(out) :: error, nc, ni
    real :: fn_cont
    external  :: fn_cont
    
    call bracketed_1d_rootfind(fn_cont, tgt_val, xl, xr, fl, fr, &
                               tol, x_soln, error,  &
                               ncalls=nc, &
                               niter=ni, &
                               absolute_convergence=1, &
                               root_convergence=1)

   end subroutine disc_abs_byroot_woutder
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function fn_cont(x)

    implicit none
    
    real :: fn_cont
    real, intent(in) :: x
    real, parameter :: pi = 3.14159265359
    real :: x_scaled
    
    x_scaled = pi * x / 2.0
    fn_cont = tan(x_scaled) - x_scaled
    
    end function fn_cont

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function dfn_cont(x)

    implicit none
    
    real :: dfn_cont
    real, intent(in) :: x
    real, parameter :: pi = 3.14159265359
    real :: x_scaled
    
    x_scaled = pi * x / 2.0
    dfn_cont = (pi/2.0) * ( 1 / (cos(x_scaled))**2 - 1.0 )
    
    end function dfn_cont

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function fn_disc(x)
  
    implicit none
    
    real :: fn_disc
    real, intent(in) :: x
    real, parameter :: disc = 0.01
    real :: x1, x2, f1, f2, slope
    real :: fn_cont
    external :: fn_cont
    
    x1 = aint(x / disc) * disc
    x2 = x1 + disc
    f1 = fn_cont(x1)
    f2 = fn_cont(x2)
    slope = ( f2 - f1 ) / ( x2 - x1 )
    fn_disc = f1 + slope * ( x - x1 )
  
  end function fn_disc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function dfn_disc(x)
  
    implicit none
    
    real :: dfn_disc
    real, intent(in) :: x
    real, parameter :: disc = 0.01
    real :: x1, x2, f1, f2, slope
    real :: fn_cont
    external :: fn_cont
    
    x1 = aint(x / disc) * disc
    x2 = x1 + disc
    f1 = fn_cont(x1)
    f2 = fn_cont(x2)
    slope = ( f2 - f1 ) / ( x2 - x1 )
    dfn_disc = slope
  
  end function dfn_disc
