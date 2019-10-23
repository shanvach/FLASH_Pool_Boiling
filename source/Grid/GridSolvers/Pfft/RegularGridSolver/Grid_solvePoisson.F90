!!****if* source/Grid/GridSolvers/Pfft/DirectSolver/Grid_solvePoisson
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
!!   Poisson solver routine.  This module implements combined 
!!   fft and matrix-based direct methods for periodic, dirichlet,
!!   and/or neumann problems on a uniform or regular grid.  
!!   Isolated problems are not supported
!!
!!   PARAMESH GRID IS NOT SUPPORTED!!!
!!
!!   Supported solver combinations are:
!!                                              fft             (2d & 3d)
!!                                              fft + tridiag   (2d & 3d)
!!                                              fft + blktri    (3d)
!!                                              fft + pdc2d     (3d)
!!                                              pdc3d           (3d)
!!
!!
!! ARGUMENTS
!!
!!  iSoln -  index to variable containing potential
!!  iSrc - index to variable containing density
!!  bcTypes - boundary types along various faces,
!!          valid values are: (although only some are implemented)
!!          GRID_PDE_BND_PERIODIC (1) (supported)
!!          GRID_PDE_BND_DIRICHLET (2) (not supported in this implementation)
!!          GRID_PDE_BND_NEUMANN (3) (supported for either the Z direction or
!!                                   for Y and Z directions, others must be periodic)
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

  use gr_pfftData, ONLY : pfft_inLen,pfft_outLen,pfft_usableProc, pfft_transformType, &
                          pfft_inLen, pfft_globalLen

  use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
                             Grid_pfftMapToInput, Grid_getGlobalIndexLimits,&
                             Grid_getBlkIndexLimits, Grid_pfftInit, 
                             Grid_pfft, Grid_pfftMapFromOutput, Grid_getBlkPtr, &
                             Grid_releaseBlkPtr, Grid_pfftFinalize

  use gr_pfftInterface, ONLY : gr_pfftPoissonDirect, gr_pfftSpecifyTransform

  use Grid_data, ONLY : gr_meshNumProcs, gr_meshMe

  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Pfft.h"  

  integer, intent(in)    :: iSoln, iSrc
  integer, intent(in)    :: bcTypes(2*MDIM)
  real, intent(in)       :: bcValues(2,2*MDIM)
  real, intent(in)       :: poisfact

  !--------------------------------------------------------------------------
  real, allocatable, dimension(:) :: inArray, tranArray, outArray
  real :: meanDensity, normScale
  integer, dimension(MDIM) :: localSize, globalSize, transformType
  integer, dimension(0:MDIM) :: baseDatType
  integer :: inSize, tranSize, solveFlag

  !=========================================================================================


  ! Define Solveflag
  ! solveflag  1 => periodic in x, Neumann conditions in y, z  (BLKTRI)
  ! solveflag  2 => periodic in x, z and Neumann in y  (TRIDIAG)
  ! solveflag  3 => periodic in x, y and z
  ! Use BC_type to define the solver, if BC arrangement not supported stop.

  if( bcTypes(1) .eq. GRID_PDE_BND_PERIODIC  .and.  &   ! in x periodic
      bcTypes(3) .eq. GRID_PDE_BND_NEUMANN   .and.  &   ! in y Neumann
      bcTypes(5) .eq. GRID_PDE_BND_NEUMANN ) then       ! in z Neumann
      
     solveflag = 1

  elseif(bcTypes(1) .eq. GRID_PDE_BND_PERIODIC .and.  & ! in x periodic
         bcTypes(3) .eq. GRID_PDE_BND_PERIODIC .and. &  ! in y periodic
         bcTypes(5) .eq. GRID_PDE_BND_NEUMANN ) then    ! in z Neumann

     solveflag = 2
 
  elseif(bcTypes(1) .eq. GRID_PDE_BND_PERIODIC .and. &   ! in x periodic
         bcTypes(3) .eq. GRID_PDE_BND_PERIODIC .and. &   ! in y periodic
         bcTypes(5) .eq. GRID_PDE_BND_PERIODIC ) then    ! in z periodic

     solveflag = 3

  else

     if (gr_meshMe .eq. 0) then
        write(*,*) 'Error: Poisson Solution Boundary Conditions Not Supported'
        write(*,*) 'bcTypes x =',bcTypes(1),bcTypes(2)
        write(*,*) 'bcTypes y =',bcTypes(3),bcTypes(4)
        write(*,*) 'bcTypes z =',bcTypes(5),bcTypes(6)
     endif
     call Driver_abortFlash("Error: Poisson Boundary Conditions Not Supported")

  endif

  ! Initialize here
  globalSize(:) = pfft_globalLen(:)
  transformType(:) = pfft_transformType(:)
  localSize(:) = pfft_inLen(:)

  !Important.  Tests that this processor should be doing work
  if(.not.pfft_usableProc) return


  inSize=pfft_inLen(IAXIS)*pfft_inLen(JAXIS)*pfft_inLen(KAXIS)
  allocate(inArray(inSize+2))
  allocate(outArray(inSize+2))


  !! Here's the real work of the fft
  ! Converts to uniform mesh (on output, inArray contains uniformly mapped density)
  call Grid_pfftMapToInput(iSrc,inArray) 

  ! Call gr_pfftPoisson Direct Forward:
  call gr_pfftPoissonDirect (PFFT_FORWARD, solveflag, inSize, localSize, globalSize, &
                           transformType, inArray, outArray)

  ! Call gr_pfftPoisson Direct Inverse: 
  call gr_pfftPoissonDirect (PFFT_INVERSE, solveflag, inSize, localSize, globalSize, & 
                           transformType, inArray, outArray)


  ! Now multiply by the poisson factor
  outArray(1:inSize) = outArray(1:inSize)*poisfact
  
  ! Map back to the non-uniform mesh
  call Grid_pfftMapFromOutput(iSoln,outArray)
  
  deallocate(inArray)
  deallocate(outArray)

  return
end subroutine Grid_solvePoisson
