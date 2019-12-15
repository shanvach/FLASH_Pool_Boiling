!!****if* source/Grid/GridSolvers/Pfft/RegularGridSolver/gr_pfftPoissonHomBcTrig
!!
!! NAME
!!
!!  gr_pfftPoissonPoissonHomBcTrig
!!
!! SYNOPSIS
!!
!!  call gr_pfftPoissonHomBcTrig(integer(IN) :: iDirection,
!!                              integer(IN) :: iSrc, 
!!                              integer(IN) :: inSize,
!!                              integer(IN) :: bcTypes(6),
!!                              integer(IN) :: bcValues(2,6),
!!                              real(IN), dimension(inSize)  :: inArray(:)
!!                              real(OUT), dimension(inSize) :: outArray(:)
!!
!! DESCRIPTION
!!
!!   Poisson solver routine.  This routine implements the simple
!!   FFT based method for periodic on a uniform grid.
!!   
!!   Isolated problems are not supported.
!!
!!
!! ARGUMENTS
!!
!!  iDirection - direction of the transform, valid values are 
!!               PFFT_FORWARD and PFFT_INVERSE
!!  iSrc       - index to variable containing density
!!  inSize     - size of inArray and outArray
!!  bcTypes    - boundary types along various faces, valid values are:
!!               GRID_PDE_BND_PERIODIC  (1) (supported)
!!               GRID_PDE_BND_DIRICHLET (2) (homogeneous Dirichlet supported)
!!               GRID_PDE_BND_NEUMANN   (3) (homogeneous Neumann supported)
!!  bcValues   - the values to boundary conditions, currently not used
!!  inArray    - single dimension array containing linearized 3D data
!!               to be transformed (compatible with pfft comm routines)
!!  outArray   - single dimension array containing transformed 
!!               linearized 3D data to be transformed
!!
!!
!!***

subroutine gr_pfftPoissonHomBcTrig (iDirection, iSrc, inSize, bcTypes, bcValues, inArray, outArray)

  use gr_pfftData, ONLY : pfft_outLen
  use Grid_interface, ONLY : Grid_pfft, GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN
  use gr_pfftInterface, ONLY : gr_pfftDerivs
  use gr_interface, ONLY:  gr_findMean

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Pfft.h"  

  integer, intent(in) :: iDirection, iSrc, inSize
  integer, intent(in) :: bcTypes(2*MDIM)
  real, intent(in)    :: bcValues(2,2*MDIM)
  real, dimension(inSize), intent(inout) :: inArray
  real, dimension(inSize), intent(out)   :: outArray
  integer :: tranSize
  real    :: meanDensity
  real, dimension(:), allocatable, save :: tranArray

 
 
  ! Forward Sweep and Solution
  if (iDirection .eq. PFFT_FORWARD) then

    ! Figure out the mean of the density and subtract
    if (ALL(bcTypes(1:2*NDIM)==GRID_PDE_BND_PERIODIC .OR. &
            bcTypes(1:2*NDIM)==GRID_PDE_BND_NEUMANN)) then
      call gr_findMean(iSrc, 2, .false., meanDensity)
      inArray(1:inSize) = inArray(1:inSize) - meanDensity
    end if

    ! Forward transform of density
    tranSize = 2 * pfft_outLen(IAXIS) * pfft_outLen(JAXIS) * pfft_outLen(KAXIS)
    allocate(tranArray(tranSize))
    call Grid_pfft(PFFT_FORWARD, inArray, tranArray)
    
    ! Calculates the transform of iSoln = GPOT
    !     which is the transform of the delSquared(u) = rho in Poisson equation
    call gr_pfftDerivs(tranArray)



  ! Inverse Transform Solution
  else

    ! Inverse transform of GPOT
    call Grid_pfft(PFFT_INVERSE,tranArray,outArray)
    deallocate(tranArray)



  end if
  
end subroutine gr_pfftPoissonPeriodic
