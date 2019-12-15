!!****if* source/Grid/localAPI/gr_pfftPoissonTrigDirect
!!
!! NAME
!!
!!  gr_pfftPoissonTrigDirect
!!
!!
!! SYNOPSIS
!!
!!  call gr_pfftPoissonTrigDirect(integer(IN)                :: iDirection,
!!                                integer(IN)                :: solveflag,
!!                                integer(IN)                :: inSize,
!!                                integer(IN)                :: localSize(MDIM),
!!                                integer(IN)                :: globalSize(MDIM),
!!                                integer(IN)                :: transformType(MDIM),
!!                                real(IN) ,dimension(inSize):: inArray(:),
!!                                real(OUT),dimension(inSize):: outArray(:))
!!
!! DESCRIPTION
!!
!!   Poisson solver routine.  This module implements the mixed fft
!!   and direct based method for periodic, neumann, and dirichlet problems
!!   on a regular grid. Mixed boundaries for a given axis are not supported.
!!
!!
!! ARGUMENTS
!!
!!   iDirection  - direction of the transform, valid values are 
!!                 PFFT_FORWARD and PFFT_INVERSE 
!!   solveflag   - Indicates the solvers to be applied.
!!                 solveflag==TRIG_DRCT      => Transform in X, Direct Solve in Y
!!                 solveflag==TRIG_TRIG_DRCT => Transform in X, Direct Solve in Y, Z
!!                 solveflag==TRI:wq
!!   inSize      - size of inArray and outArray
!!   localSize   - the local bounds (on myPe) of the original 3D data to be transformed
!!   globalSize  - global size of the 3D data to be transformed
!!   transformType - The type if transform to be applied along each of the dimensions
!!   inArray       - single dimension array containing linearized data to be transformed
!!   outArray      - single dimension array containing transformed linearized data
!!  
!!***

subroutine gr_pfftPoissonTrigDirect (iDirection, solveflag, inSize, localSize, globalSize, transformType, inArray, outArray)
  
#include "constants.h"  
  
  implicit none
  
  integer, intent(in)    :: iDirection, solveflag, inSize  
  integer, dimension(MDIM),intent(IN) :: localSize,globalSize,transformType
  real,  dimension(inSize),intent(IN) :: inArray
  real,  dimension(inSize),intent(OUT) :: outArray
  outArray = 0.0
  return
end subroutine gr_pfftPoissonTrigDirect


