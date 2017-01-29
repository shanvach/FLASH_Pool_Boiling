!*******************************************************************************

!  Routine:     pois_solve_2d

!  Description: Solve Poisson's equation on a 2D rectangular domain using
!               Fourier, sine, or cosine transforms.  The type of transform
!               employed depends upon the requested boundary conditions.
!               Transforms are performed using Takuya Ooura's FFT package
!               (http://momonga.t.u-tokyo.ac.jp/~ooura/fft.html).

!  Parameters:  soln(:,:)     On input, source array, sampled at zone centers;
!                             on output, array to receive solution, sampled at
!                             zone centers
!               nx, ny        Dimensions of source and solution arrays
!               dx, dy        Zone widths in x, y, and z
!               ibnd          Type of boundary conditions to assume:
!                             0 = periodic, 1 = Dirichlet, 2 = Neumann

module Module2DCGPoissonSolver

contains

subroutine pois_solve_2d (soln, nx, ny, dx, dy, ibnd, &
                          Gx, Gy, wk1, wk2, ip, norm, nwk1, nwk2, nip, G2D)

!===============================================================================

real, intent(in)       :: dx, dy, norm
integer, intent(in)    :: nx, ny, ibnd, nwk1, nwk2, nip
real, intent(inout)    :: soln(0:nx-1,0:ny-1)
real, intent(inout)    :: wk1(0:nwk1), wk2(0:nwk2), Gx(0:nx), Gy(0:ny)
integer, intent(inout) :: ip(0:nip)
real, intent(inout)    :: G2D(0:nx,0:ny)

integer         :: i, j
real            :: G

!===============================================================================

! Solve the Poisson equation by transforming the source, applying the Green's
! function, and performing the inverse transform.

select case (ibnd)

!===============================================================================

! Periodic boundary conditions:  use real discrete Fourier transform

  case (0)

    call rdft2d (nx, nx, ny, 1, soln, wk1, ip, wk2)

    do j = 1, ny-1
      do i = 1, nx/2-1
!        G = GreenFctn2D(Gx(i), Gy(j), ibnd)
        G = G2D(i,j)
        soln(2*i,j)   = soln(2*i,j) * G
        soln(2*i+1,j) = soln(2*i+1,j) * G
      enddo
    enddo

    do i = 1, nx/2-1
!      G = GreenFctn2D(Gx(i), Gy(0), ibnd)
      G = G2D(i,0)
      soln(2*i,0)   = soln(2*i,0) * G
      soln(2*i+1,0) = soln(2*i+1,0) * G
    enddo

    do j = 1, ny/2-1
!      G = GreenFctn2D(Gx(0), Gy(j), ibnd)
      G = G2D(0,j)
      soln(0,j) = soln(0,j) * G
      soln(1,j) = soln(1,j) * G
!      G = GreenFctn2D(Gx(nx/2), Gy(j), ibnd)
      G = G2D(nx/2,j)
      soln(1,ny-j) = soln(1,ny-j) * G
      soln(0,ny-j) = soln(0,ny-j) * G
    enddo

!    G = GreenFctn2D(Gx(0), Gy(0), ibnd)
    G = G2D(0,0)
    soln(0,0) = soln(0,0) * G
!    G = GreenFctn2D(Gx(nx/2), Gy(0), ibnd)
    G = G2D(nx/2,0)
    soln(1,0) = soln(1,0) * G
!    G = GreenFctn2D(Gx(0), Gy(ny/2), ibnd)
    G = G2D(0,ny/2)
    soln(0,ny/2) = soln(0,ny/2) * G
!    G = GreenFctn2D(Gx(nx/2), Gy(ny/2), ibnd)
    G = G2D(nx/2,ny/2)
    soln(1,ny/2) = soln(1,ny/2) * G

    call rdft2d (nx, nx, ny, -1, soln, wk1, ip, wk2)

    do j = 0, ny-1
      do i = 0, nx-1
        soln(i,j) = soln(i,j) * norm
      enddo
    enddo

!===============================================================================

! Dirichlet boundary conditions:  use discrete sine transform

  case (1)

    call ddst2d (nx, nx, ny, -1, soln, wk1, ip, wk2)

    do j = 1, ny
      do i = 1, nx
!        G = GreenFctn2D(Gx(i), Gy(j), ibnd)
        G = G2D(i,j)
        soln(mod(i,nx),mod(j,ny)) = soln(mod(i,nx),mod(j,ny)) * G
      enddo
    enddo

    do i = 0, nx-1
      soln(i,0) = soln(i,0) * 0.5
    enddo

    do j = 0, ny-1
      soln(0,j) = soln(0,j) * 0.5
    enddo

    call ddst2d (nx, nx, ny, 1, soln, wk1, ip, wk2)

    do j = 0, ny-1
      do i = 0, nx-1
        soln(i,j) = soln(i,j) * norm
      enddo
    enddo

!===============================================================================

! Neumann boundary conditions:  use discrete cosine transform

  case (2)

    call ddct2d (nx, nx, ny, -1, soln, wk1, ip, wk2)

    do j = 0, ny-1
      do i = 0, nx-1
!        G = GreenFctn2D(Gx(i), Gy(j), ibnd)
        G = G2D(i,j)
        soln(i,j) = soln(i,j) * G
      enddo
    enddo

    do i = 0, nx-1
      soln(i,0) = soln(i,0) * 0.5
    enddo

    do j = 0, ny-1
      soln(0,j) = soln(0,j) * 0.5
    enddo

    call ddct2d (nx, nx, ny, 1, soln, wk1, ip, wk2)

    do j = 0, ny-1
      do i = 0, nx-1
        soln(i,j) = soln(i,j) * norm
      enddo
    enddo

!===============================================================================

end select

!===============================================================================

return
end subroutine pois_solve_2d

!*******************************************************************************

!  Routine:     GreenFctn2D

!  Description: Evaluate the Green's function for the Poisson equation at a
!               specified k-value.  x and y components of the denominator
!               sum are accepted as arguments.


function GreenFctn2D (Gx, Gy, ibnd)

!===============================================================================

implicit none

real, intent(in)    :: Gx, Gy
integer, intent(in) :: ibnd
real                :: GreenFctn2D

real                :: Ginv

!===============================================================================

Ginv = Gx + Gy

if (abs(Ginv) > 1.E-99) then
  GreenFctn2D = 1./Ginv
else
  GreenFctn2D = 0.
endif

!===============================================================================

return
end function GreenFctn2D

end module Module2DCGPoissonSolver
