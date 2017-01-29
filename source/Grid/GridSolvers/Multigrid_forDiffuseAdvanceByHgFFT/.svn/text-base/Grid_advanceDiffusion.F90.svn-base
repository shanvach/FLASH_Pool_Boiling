!!****if* source/Grid/GridSolvers/Multigrid_forDiffuseAdvanceByHgFFT/Grid_advanceDiffusion
!!
!! NAME
!!  Grid_advanceDiffusion
!!
!! SYNOPSIS
!!
!!  call Grid_advanceDiffusion(integer(IN) :: iSoln,
!!                             integer(IN) :: iSrc,
!!                             integer(IN) :: iCond,
!!                             integer(IN) :: bcTypes(6),
!!                             real(IN)    :: bcValues(2,6),
!!                             real(IN)    :: dt,
!!                             real(IN)    :: chi,
!!                             real(IN)    :: scaleFact,
!!                             logical(IN) :: solnIsDelta)
!!
!! DESCRIPTION
!!
!!  This routine advances a diffusive operator (think head conduction) in time
!!  for assorted boundary types using the Huang-Greengard multilevel solver framework.
!!  
!! ARGUMENTS
!!
!!  iSoln    - the index for the solution variable (potential when used for self-gravity)
!!  iSrc     - the index of the source variable (density when used for self-gravity)
!!  iCond    - the index of the diffusion coefficient (if used, otherwise -1)
!!  bcTypes  - the boundary condition type; only the first entry is used.
!!                    the following conditions are presently supported:
!!                    MG_BND_PERIODIC
!!                    MG_BND_ISOLATED
!!                    MG_BND_DIRICHLET
!!                    MG_BND_NEUMANN could most likely be easily implemented
!!  bcValues - an unused argument used to keep the interface standard
!!  dt             - the time step
!!  chi            - factor in the heat diffusion equation, here assumed constant (!)
!!                       d/dt T  =  chi * Laplace T
!!  scaleFact      - the scaling factor for the eventual solution
!!  solnIsDelta    - Is the solution only a delta that the caller has to apply to the
!!                   temperature, rather than temperature itself?
!!
!! RESULT
!!  
!!  The variable at iSoln contains the solution to L(u) = rho for the
!!  given boundary value setup
!!
!! PARAMETERS
!!  ISLS_VAR -- residual grid variable
!!  ICOR_VAR -- correction grid variable
!!
!! NOTES
!!
!!  It is currently assumed in this implementation
!!  DEV: But not checked! - KW
!!  that boundary condition types at all sides of the domain are the same. 
!!  What this means to the user  is that only the first value in bcTypes is
!!  checked here, and may be assumed to give the type of boundary condition
!!  for all 2*NDIM directions.
!!  Presently, the FFT package used only supports the same transform done in
!!  all directions.  Using something different, such as PFFT as the top level,
!!  should remedy this.  
!!  The all-sides-the-same boundary values could be for the sake of the
!!  top-level FFT solve (but the dummy argument is not actually used in this
!!  version).
!!
!! SEE ALSO
!! 
!!  gr_hgSolve
!!  
!!
!!***

subroutine Grid_advanceDiffusion (iSoln, iSrc, iCond, bcTypes, bcValues, dt, chi, scaleFact, &
     solnIsDelta,pass)

  use gr_hgData, ONLY : gr_hgMeshRefineMax
  use Grid_data, ONLY : gr_meshMe

  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_hgInterface, ONLY: gr_hgSolve

  implicit none

#include "Multigrid.h"
#include "Flash.h"
#include "constants.h"

  integer, intent(IN)               :: iSoln, iSrc
  integer, intent(IN)               :: iCond
  integer, dimension(6), intent(IN) :: bcTypes
  real, dimension(2,6), intent(IN)  :: bcValues
  real, intent(IN)                  :: dt, chi, scaleFact
  logical, intent(IN)               :: solnIsDelta
  integer, OPTIONAL, intent(IN)     :: pass
  real                              :: norm

  integer               :: i

  external              gr_hgDiffuseSolveBlock

  !=======================================================================


  ! Initializations performed in the stand-alone Diffuse unit to be.

  ! Call the Huang-Greengard multigrid solver framework for different types of boundary
  ! conditions.

  select case (bcTypes(1))

     !---------------------------------------------------------------

  case (MG_BND_ISOLATED)

     
     call Driver_abortFlash("[Grid_advanceDiffusion]  MG_BND_ISOLATED shall only be used for the Poisson solver!")
     
     !-------------------------------------------------------------------------------
     
  case (MG_BND_PERIODIC)
     
     call Timers_start("Multigrid_solve")

     call gr_hgSolve(iSrc, iSoln, ISLS_VAR, ICOR_VAR, gr_hgDiffuseSolveBlock, &
          bcTypes, scaleFact, dt, chi)

     call Timers_stop("Multigrid_solve")

     !----------------------------------------------------------------

  case (MG_BND_DIRICHLET)

     call Timers_start("Multigrid_solve")

     call gr_hgSolve(iSrc, iSoln, ISLS_VAR, ICOR_VAR, gr_hgDiffuseSolveBlock, &
          bcTypes, scaleFact, dt, chi)

     call Timers_stop("Multigrid_solve")

     !----------------------------------------------------------------------

  case (MG_BND_NEUMANN)

     call Timers_start("Multigrid_solve")

     call gr_hgSolve(iSrc, iSoln, ISLS_VAR, ICOR_VAR, gr_hgDiffuseSolveBlock, &
          bcTypes, scaleFact, dt, chi)

     call Timers_stop("Multigrid_solve")

     !----------------------------------------------------------------------

  case default

     call Driver_abortFlash("[Grid_advanceDiffusion]  invalid boundary condition type!")

     !-------------------------------------------------------------------------------

  end select
  !===============================================================================



 return
end subroutine Grid_advanceDiffusion
