!*******************************************************************************

! Routine:      poisson_hg_solve_block

! Description:  Single-block Poisson solver for use with the Huang-Greengard
!               multilevel adaptive solver.

!               This version uses fast Fourier transforms (FFTs) to solve
!               the Poisson equation.

! Parameters:   soln(:,:,:)   On input, source array, sampled at zone centers;
!                             on output, array to receive solution, sampled at
!                             zone centers
!               nx, ny, nz    Dimensions of source and solution arrays
!               dx, dy, dz    Zone widths in x, y, and z
!               bnd_type      Type of boundary conditions to assume:
!                             0 = periodic, 1 = Dirichlet, 2 = Neumann


subroutine poisson_hg_solve_block (soln, nx, ny, nz, dx, dy, dz, bnd_type,&
                                   level)

!===============================================================================

use runtime_parameters
use dBase, ONLY: ndim

use Module1DCGPoissonSolver
use Module2DCGPoissonSolver
use Module3DCGPoissonSolver

implicit none

integer, intent(in) :: nx, ny, nz, bnd_type, level
real, intent(in)    :: dx, dy, dz
real, intent(inout) :: soln(nx,ny,nz)

logical, save              :: first_call = .true.
integer, save              :: lrefine_max, nwk1, nwk2, nip
integer, allocatable, save :: iinit(:), ip(:,:)
real, allocatable, save    :: wk1(:,:), wk2(:,:), Gx(:,:), Gy(:,:), Gz(:,:)
real, allocatable, save    :: norm(:), G2D(:,:,:), G3D(:,:,:,:)
real, parameter            :: pi = 3.1415926535898

real                       :: nxinv, nyinv, nzinv, dxinv, dyinv, dzinv
integer                    :: m, i, j, k

!===============================================================================

if (first_call) then
  call get_parm_from_context("lrefine_max", lrefine_max)
  allocate(iinit(lrefine_max))
  iinit(:) = 1
  allocate(norm(lrefine_max))
  allocate(Gx(0:nx,lrefine_max))
  allocate(Gy(0:ny,lrefine_max))
  allocate(Gz(0:nz,lrefine_max))
  if (ndim == 2) then
    allocate(G2D(0:nx,0:ny,lrefine_max))
  elseif (ndim == 3) then
    allocate(G3D(0:nx,0:ny,0:nz,lrefine_max))
  endif
  if (ndim == 1) then
    nwk1 = nx*5/4-1
    nwk2 = 0
    nip  = 1+int(sqrt(float(nx/2)))
  else
!    select case (bnd_type)
!      case (0) ! periodic
        nwk1 = 8*max(ny, nz)-1
!        nwk2 = max(nx/4, ny/2, nz/2)+nx/4-1
        nip  = 1+int(sqrt(float(max(nx/2, ny, nz))))
!      case (1) ! Dirichlet
!        nwk1 = 4*max(ny, nz)-1
        nwk2 = max(nx*3/2, ny*3/2, nz*3/2)-1
!        nip  = 1+int(sqrt(float(max(nx/2, ny/2, nz/2))))
!      case (2) ! Neumann
!        nwk1 = 4*max(ny, nz)-1
!        nwk2 = max(nx*3/2, ny*3/2, nz*3/2)-1
!        nip  = 1+int(sqrt(float(max(nx/2, ny/2, nz/2))))
!    end select
  endif
  allocate (wk1(0:nwk1,lrefine_max))
  allocate (wk2(0:nwk2,lrefine_max))
  allocate (ip(0:nip,lrefine_max))
  first_call = .false.
endif

if (iinit(level) == 1) then
  dxinv = 2./dx**2
  dyinv = 2./dy**2
  dzinv = 2./dz**2
  nxinv = pi/nx
  nyinv = pi/ny
  nzinv = pi/nz
  Gx(:,level) = 0.
  Gy(:,level) = 0.
  Gz(:,level) = 0.
  select case (bnd_type)
    case (0) ! periodic
      do i = 0, nx-1
        Gx(i,level) = dxinv * (cos(2.*i*nxinv) - 1.)
      enddo
      do j = 0, ny-1
        Gy(j,level) = dyinv * (cos(2.*j*nyinv) - 1.)
      enddo
      do k = 0, nz-1
        Gz(k,level) = dzinv * (cos(2.*k*nzinv) - 1.)
      enddo
      norm(level) = 2./(float(nx)*float(ny)*float(nz))
    case (1) ! Dirichlet
      do i = 1, nx
        Gx(i,level) = dxinv * (cos(i*nxinv) - 1.)
      enddo
      do j = 1, ny
        Gy(j,level) = dyinv * (cos(j*nyinv) - 1.)
      enddo
      do k = 1, nz
        Gz(k,level) = dzinv * (cos(k*nzinv) - 1.)
      enddo
      norm(level) = 2.**ndim / (float(nx)*float(ny)*float(nz))
    case (2) ! Neumann
      do i = 0, nx-1
        Gx(i,level) = dxinv * (cos(i*nxinv) - 1.)
      enddo
      do j = 0, ny-1
        Gy(j,level) = dyinv * (cos(j*nyinv) - 1.)
      enddo
      do k = 0, nz-1
        Gz(k,level) = dzinv * (cos(k*nzinv) - 1.)
      enddo
      norm(level) = 2.**ndim / (float(nx)*float(ny)*float(nz))
  end select
  if (ndim == 2) then
    do j = 0, ny
      do i = 0, nx
        if (abs(Gx(i,level)+Gy(j,level)) > 1.E-99) then
          G2D(i,j,level) = 1./(Gx(i,level)+Gy(j,level))
        else
          G2D(i,j,level) = 0.
        endif
      enddo
    enddo
  elseif (ndim == 3) then
    do k = 0, nz
      do j = 0, ny
        do i = 0, nx
          if (abs(Gx(i,level)+Gy(j,level)+Gz(k,level)) > 1.E-99) then
            G3D(i,j,k,level) = 1./(Gx(i,level)+Gy(j,level)+Gz(k,level))
          else
            G3D(i,j,k,level) = 0.
          endif
        enddo
      enddo
    enddo
  endif
  ip(0,level) = 0
  iinit(level) = 0
endif

if (ndim == 1) then
  call pois_solve_1d (soln, nx, dx, bnd_type, &
                      Gx(0,level), wk1(0,level), &
                      ip(0,level), norm(level), nwk1, nip)
elseif (ndim == 2) then
  call pois_solve_2d (soln, nx, ny, dx, dy, bnd_type, &
                      Gx(0,level), Gy(0,level), wk1(0,level), &
                      wk2(0,level), ip(0,level), norm(level), &
                      nwk1, nwk2, nip, G2D(0,0,level))
else
  call pois_solve_3d (soln, nx, ny, nz, dx, dy, dz, bnd_type, &
                      Gx(0,level), Gy(0,level), Gz(0,level), &
                      wk1(0,level), wk2(0,level), ip(0,level), norm(level), &
                      nwk1, nwk2, nip, G3D(0,0,0,level))
endif

!===============================================================================

return
end subroutine poisson_hg_solve_block
