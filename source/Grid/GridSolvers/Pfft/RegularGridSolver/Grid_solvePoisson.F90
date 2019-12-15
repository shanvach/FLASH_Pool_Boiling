!!****if* source/Grid/GridSolvers/Pfft/RegularGridSolver/Grid_solvePoisson
!!
!! NAME
!!
!!  Grid_solvePoisson
!!
!! SYNOPSIS
!!
!!  call Grid_solvePoisson(integer(IN) :: iSoln,
!!                         integer(IN) :: iSrc, 
!!                         integer(IN) :: bcTypes(6),
!!                         real(IN)    :: bcValues(2,6),
!!                         real(INOUT) :: poisfact)
!!
!! DESCRIPTION
!!
!!   Poisson solver routine.  This implementation provides the
!!   FFT based method for periodic, homogeneous Dirichlet, and
!!   homogeneous Neumann problems on a regular grid.
!!   Mixed boundary problems, in which different boundary types
!!   apply at different boundary faces, are also supported.
!!   
!!   Isolated problems are not supported.
!!
!!
!! ARGUMENTS
!!
!!  iSoln -  index to variable containing potential
!!  iSrc - index to variable containing density
!!  bcTypes - boundary types along various faces,
!!          valid values are: (although only some are implemented)
!!          GRID_PDE_BND_PERIODIC (1) (supported)
!!          GRID_PDE_BND_DIRICHLET (2) (homogeneous Dirichlet supported)
!!          GRID_PDE_BND_NEUMANN (3) (homogeneous Neumann supported)
!!          GRID_PDE_BND_ISOLATED (0) (not supported in this implementation)
!!  bcValues - the values to boundary conditions, currently not used
!!  poisfact -  factor to be used in calculation
!!
!! NOTES
!!
!!  The symbols listed above for bcTypes are declared as FORTRAN PARAMETERS in
!!  the module Grid_interfaces.  Code using this interface should refer to that
!!  module with a USE statement, like this:
!!
!!    use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
!!       GRID_PDE_BND_DIRICHLET, &
!!       Grid_solvePoisson
!!
!!***

subroutine Grid_solvePoisson (iSoln, iSrc, bcTypes, bcValues, poisfact)

  use gr_pfftData, ONLY : pfft_inLen, pfft_usableProc, pfft_globalLen, pfft_transformType, pfft_solver
  use Grid_interface, ONLY : Grid_pfftMapToInput, Grid_pfftMapFromOutput, GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN
  use gr_pfftInterface, ONLY : gr_pfftPoissonHomBcTrig, gr_pfftPoissonTrigDirect, gr_pfftPoissonDirect


  implicit none

#include "constants.h"
#include "Flash.h"
#include "Pfft.h"  

  integer, intent(in)    :: iSoln, iSrc
  integer, intent(in)    :: bcTypes(2*MDIM)
  real, intent(in)       :: bcValues(2,2*MDIM)
  real, intent(inout)    :: poisfact 
  real, allocatable, dimension(:) :: inArray, outArray
  integer, dimension(MDIM) :: localSize, globalSize, transformType
  integer :: inSize
  real    :: meanDensity = 0.

  !Important.  Tests that this processor should be doing work
  if(.not.pfft_usableProc) return   

  inSize=pfft_inLen(IAXIS)*pfft_inLen(JAXIS)*pfft_inLen(KAXIS)
  allocate(inArray(inSize+2))
  allocate(outArray(inSize+2))

  inArray(:) = 0.
  outArray(:) = 0.

  globalSize(:) = pfft_globalLen(:)
  transformType(:) = pfft_transformType(:)
  localSize(:) = pfft_inLen(:)

  ! map to uniform mesh 
  call Grid_pfftMapToInput(iSrc,inArray) 

  ! use the appropriate poisson solver
  select case (pfft_solver)

  ! Homogenious boundary conditions w/o stretching in {x, y}
  case (TRIG_TRIG, TRIG_TRIG_TRIG)

    ! Call gr_pfftPoisson Periodic Forward:
    call gr_pfftPoissonHomBcTrig (PFFT_FORWARD, iSrc, inSize, bcTypes, bcValues, inArray, outArray)
  
    ! Call gr_pfftPoisson Periodic Inverse:
    call gr_pfftPoissonHomBcTrig (PFFT_INVERSE, iSrc, inSize, bcTypes, bcValues, inArray, outArray)
  
    ! Now multiply by the poisson factor
    outArray(1:inSize) = outArray(1:inSize)*poisfact

  ! Homogenious boundary conditions w/o x-axis stretching in {x, y}
  case (TRIG_DRCT, TRIG_TRIG_DRCT, TRIG_DRCT_DRCT)

    ! Call gr_pfftPoisson Periodic Forward:
    call gr_pfftPoissonTrigDirect (PFFT_FORWARD, 1, inSize, localSize, globalSize, transformType, inArray, outArray)
  
    ! Call gr_pfftPoisson Periodic Inverse:
    call gr_pfftPoissonTrigDirect (PFFT_INVERSE, 1, inSize, localSize, globalSize, transformType, inArray, outArray)
  
    ! Now multiply by the poisson factor
    outArray(1:inSize) = outArray(1:inSize)*poisfact

  ! Homogenious boundary conditions w/ stretching in {x, y}
  case (DRCT_DRCT, DRCT_DRCT_DRCT)

    ! Call gr_pfftPoisson Periodic Forward:
    call gr_pfftPoissonDirect (PFFT_FORWARD, 1, inSize, localSize, globalSize, transformType, inArray, outArray)
  
    ! Call gr_pfftPoisson Periodic Inverse:
    call gr_pfftPoissonDirect (PFFT_INVERSE, 1, inSize, localSize, globalSize, transformType, inArray, outArray)
  
    ! Now multiply by the poisson factor
    outArray(1:inSize) = outArray(1:inSize)*poisfact


  case default
    call Driver_abortFlash("Specification of Poisson Solver options not supported!")
 
  end select

  ! Map back to the non-uniform mesh
  call Grid_pfftMapFromOutput(iSoln,outArray)

  deallocate(inArray)
  deallocate(outArray)

  return
end subroutine Grid_solvePoisson
