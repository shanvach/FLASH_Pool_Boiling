!!****if* source/Grid/localAPI/gr_pfftPoissonHomBcTrig
!!
!! NAME
!!
!!  gr_pfftPoissonHomBcTrig
!!
!!
!! SYNOPSIS
!!
!!  call gr_pfftPoissonHomBcTrig(integer(IN) :: iDirection,
!!                               integer(IN) :: iSrc, 
!!                               integer(IN) :: inSize,
!!                               integer(IN) :: bcTypes(6),
!!                               integer(IN) :: bcValues(2,6),
!!                               real(IN), dimension(inSize)  :: inArray(:)
!!                               real(OUT), dimension(inSize) :: outArray(:)
!!
!! DESCRIPTION
!!
!!   Poisson solver routine.  This routine implements the simple
!!   FFT based method for periodic, diriclet, and neumann on a uniform grid.
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

  implicit none

#include "constants.h"

  integer, intent(in)    :: iDirection, iSrc, inSize
  integer, intent(in) :: bcTypes(2*MDIM)
  real, intent(in)    :: bcValues(2,2*MDIM)
  real, dimension(inSize), intent(inout) :: inArray
  real, dimension(inSize), intent(out)   :: outArray

  outArray = 0.0
  return

end subroutine gr_pfftPoissonHomBcTrig


