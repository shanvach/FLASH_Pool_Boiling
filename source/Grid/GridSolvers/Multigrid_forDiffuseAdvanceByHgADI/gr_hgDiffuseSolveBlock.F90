!!****if* source/Grid/GridSolvers/Multigrid_forDiffuseAdvanceByHgADI/gr_hgDiffuseSolveBlock
!!
!! NAME
!!  gr_hgDiffuseSolveBlock
!!
!! SYNOPSIS
!!  gr_hgDiffuseSolveBlock(real, intent(inout) :: soln(:,:,:),
!!                         integer, intent(in) :: nx,
!!                         integer, intent(in) :: ny,
!!                         integer, intent(in) :: nz,
!!                         real, intent(in)    :: dx,
!!                         real, intent(in)    :: dy,
!!                         real, intent(in)    :: dz,
!!                         integer, intent(in) :: bnd_type,
!!                         integer, intent(in) :: level,
!!                         real(IN)            :: dt,
!!                         real(IN)            :: chi,
!!                         real(IN)            :: theta,
!!                         real(IN)            :: oldSrc,
!!                         integer(IN)         :: sol_ile,sol_iue,sol_jle,sol_jue,sol_kle,sol_kue,
!!                         integer(IN)         :: src_ile,src_iue,src_jle,src_jue,src_kle,src_kue)
!!
!! DESCRIPTION
!!  
!!  Single block Diffuse Solver using a Crank Nicholson scheme.
!!
!! ARGUMENTS
!!
!!  soln       - the input and output arrays, sampled at zone centers
!!  nx, ny, nz - the dimensions of the array
!!  dx, dy, dz - the zone widths
!!  bnd_type   - MG_BND_PERIODIC, MG_BND_DIRICHLET, or MG_BND_NEUMANN
!!  level      - the current level to solve at (precalculation of dimensions)
!!  dt         - time step, to be passed down
!!  chi        - a factor, to be passed down 
!!  theta      - a factor to switch between Implicit/Explicit Scheme
!!  oldSrc     - the old source variable 
!!  sol_ile,sol_iue,sol_jle,sol_jue,sol_kle,sol_kue - lower and upper bounds of solution array
!!  src_ile,src_iue,src_jle,src_jue,src_kle,src_kue - lower and upper bounds of source array
!!
!! NOTES
!!
!! This routine is called as SolveBlock from gr_hgSolveLevel with arguments
!!    (soln,NXB,NYB,NZB,dx,dy,dz,bnd_type,level,dt, chi, theta, oldSrc,sol_ile, &
!!     sol_iue,sol_jle,sol_jue,sol_kle,sol_kue,src_ile,src_iue,src_jle,src_jue,src_kle,src_kue)
!!
!!***

subroutine gr_hgDiffuseSolveBlock (soln, nx, ny, nz, dx, dy, dz, bnd_type,&
                                   level, &
                                   dt, chi, theta, oldSrc,&
                                   sol_ile,sol_iue,sol_jle,sol_jue,sol_kle,sol_kue,&
                                   src_ile,src_iue,src_jle,src_jue,src_kle,src_kue)

!=========================================================

!use runtime_parameters

#include "Multigrid.h"
#include "Flash.h"

#ifdef FLASH_GRID_PARAMESH2
use tree, ONLY: lrefine_max,ndim
#else
use paramesh_dimensions, ONLY: ndim
use tree, ONLY: lrefine_max
#endif

  implicit none

  integer, intent(in) :: NX, NY, NZ, bnd_type, level

  real, intent(in)    :: dx, dy, dz
  real, intent(in)    :: dt, chi, theta

  integer, intent(in) :: sol_ile,sol_iue,sol_jle,sol_jue,sol_kle,sol_kue
  integer, intent(in) :: src_ile,src_iue,src_jle,src_jue,src_kle,src_kue

  real, intent(inout) :: soln(sol_ile:sol_iue,sol_jle:sol_jue,sol_kle:sol_kue)
  real, intent(in)    :: oldSrc(src_ile:src_iue,src_jle:src_jue,src_kle:src_kue)

  integer                    :: i, j, k

  real                       :: Cond_L,Cond_R, R
  real                       :: AA(nx), DD(nx), BB(nx)

  real :: L2Err = 1.0, TOL = 1.0E-12, OldTemp
  real :: LD, MD, UD
  !===============================================================================

  ! Relaxation based approach.

  L2Err = 1.0
   
  do while (sqrt(L2Err) .ge. TOL)

     L2Err = 0.0

     do k = 1, NZ
        do j = 1, NY
           do i = 1, NX

              OldTemp = soln(i,j,k)

              Cond_R = chi ! TO DO : Use COND_VAR
              Cond_L = chi ! TO DO : Use COND_VAR

              LD = -theta*Cond_L*(dt/(dx**2))
              UD = -theta*Cond_R*(dt/(dx**2))

              MD = 1.0 + theta*(Cond_R + Cond_L)*(dt/(dx**2))

              soln(i,j,k) = oldSrc(i,j,k) - (LD*soln(i-1,j,k) + UD*soln(i+1,j,k))               

              if (NDIM .ge. 2) then

                  Cond_R = chi ! TO DO : Use COND_VAR
                  Cond_L = chi ! TO DO : Use COND_VAR

                  MD = MD + theta*(Cond_R + Cond_L)*(dt/(dy**2))

                  LD = -theta*Cond_L*(dt/(dy**2))
                  UD = -theta*Cond_R*(dt/(dy**2))

                  soln(i,j,k) = soln(i,j,k) - (LD*soln(i,j-1,k) + UD*soln(i,j+1,k))

             endif

             if (NDIM .ge. 3) then

                 Cond_R = chi ! TO DO : Use COND_VAR
                 Cond_L = chi ! TO DO : Use COND_VAR

                 MD = MD + theta*(Cond_R + Cond_L)*(dt/(dz**2))

                 LD = -theta*Cond_L*(dt/(dz**2))
                 UD = -theta*Cond_R*(dt/(dz**2))

                  soln(i,j,k) = soln(i,j,k) - (LD*soln(i,j,k-1) + UD*soln(i,j,k+1))

             endif

             soln(i,j,k) = soln(i,j,k) / MD

             L2Err = L2Err + ((OldTemp-soln(i,j,k))/OldTemp)**2

           end do   
        end do
     end do 

  end do


#ifdef ELIMINATION
#endif
  !===============================================================================

  return
end subroutine gr_hgDiffuseSolveBlock
