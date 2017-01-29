!!****ih* source/flashUtilities/1dRootfind/ut_bracketed_1d_rootfindinterface
!!
!! This is the public interface file for the 1d bracketed rootfinder.
!!
!!***
module ut_bracketed_1d_rootfindinterface

  implicit none
  
  interface

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

    end subroutine bracketed_1d_rootfind

  end interface

end module ut_bracketed_1d_rootfindinterface

