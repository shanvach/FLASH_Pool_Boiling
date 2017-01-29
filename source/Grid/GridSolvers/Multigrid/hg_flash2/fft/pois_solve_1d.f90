!*******************************************************************************

!  Routine:     pois_solve_1d

!  Description: Solve Poisson's equation on a 1D rectangular domain using
!               Fourier, sine, or cosine transforms.  The type of transform
!               employed depends upon the requested boundary conditions.
!               Transforms are performed using Takuya Ooura's FFT package
!               (http://momonga.t.u-tokyo.ac.jp/~ooura/fft.html).

!  Parameters:  soln(:)       On input, source array, sampled at zone centers;
!                             on output, array to receive solution, sampled at
!                             zone centers
!               nx, ny        Dimensions of source and solution arrays
!               dx, dy        Zone widths in x, y, and z
!               ibnd          Type of boundary conditions to assume:
!                             0 = periodic, 1 = Dirichlet, 2 = Neumann

module Module1DCGPoissonSolver

contains

subroutine pois_solve_1d (soln, nx, dx, ibnd, &
                          Gx, wk1, ip, norm, nwk1, nip)

!===============================================================================

real, intent(in)       :: dx, norm
integer, intent(in)    :: nx, ibnd, nwk1, nip
real, intent(inout)    :: soln(0:nx-1)
real, intent(inout)    :: wk1(0:nwk1), Gx(0:nx)
integer, intent(inout) :: ip(0:nip)

integer         :: i
real            :: G

!===============================================================================

! Solve the Poisson equation by transforming the source, applying the Green's
! function, and performing the inverse transform.

select case (ibnd)

!===============================================================================

! Periodic boundary conditions:  use real discrete Fourier transform

  case (0)

    call rdft (nx, 1, soln, ip, wk1)

    do i = 0, nx/2-1
      G = GreenFctn1D(Gx(i), ibnd)
      soln(2*i) = soln(2*i) * G
    enddo

    do i = 1, nx/2-1
      G = GreenFctn1D(Gx(i), ibnd)
      soln(2*i+1) = soln(2*i+1) * G
    enddo

    G = GreenFctn1D(Gx(nx/2), ibnd)
    soln(1) = soln(1) * G

    call rdft (nx, -1, soln, ip, wk1)

    do i = 0, nx-1
      soln(i) = soln(i) * norm
    enddo

!===============================================================================

! Dirichlet boundary conditions:  use discrete sine transform

  case (1)

    call ddst (nx, -1, soln, ip, wk1)

    do i = 1, nx
      G = GreenFctn1D(Gx(i), ibnd)
      soln(mod(i,nx)) = soln(mod(i,nx)) * G
    enddo

    soln(0) = soln(0) * 0.5

    call ddst (nx, 1, soln, ip, wk1)

    do i = 0, nx-1
      soln(i) = soln(i) * norm
    enddo

!===============================================================================

! Neumann boundary conditions:  use discrete cosine transform

  case (2)

    call ddct (nx, -1, soln, ip, wk1)

    do i = 0, nx-1
      G = GreenFctn1D(Gx(i), ibnd)
      soln(i) = soln(i) * G
    enddo

    soln(0) = soln(0) * 0.5

    call ddct (nx, 1, soln, ip, wk1)

    do i = 0, nx-1
      soln(i) = soln(i) * norm
    enddo

!===============================================================================

end select

!===============================================================================

return
end subroutine pois_solve_1d

!*******************************************************************************

!  Routine:     GreenFctn1D

!  Description: Evaluate the Green's function for the Poisson equation at a
!               specified k-value.  The denominator is accepted as an argument.


function GreenFctn1D (Gx, ibnd)

!===============================================================================

implicit none

real, intent(in)    :: Gx
integer, intent(in) :: ibnd
real                :: GreenFctn1D

real                :: Ginv

!===============================================================================

Ginv = Gx

if (abs(Ginv) > 1.E-99) then
  GreenFctn1D = 1./Ginv
else
  GreenFctn1D = 0.
endif

!===============================================================================

return
end function GreenFctn1D

end module Module1DCGPoissonSolver
