!!****if* source/Simulation/SimulationMain/unitTest/Pfft_TransposeTest/gr_pfftSpecifyTransform
!!
!! NAME
!!
!! gr_pfftSpecifyTransform
!!
!! SYNOPSIS
!!  
!! gr_pfftSpecifyTransform(integer(OUT) :: transformType(MDIM))
!! 
!! DESCRIPTION
!!
!! Initialises the transform type for the particular solver.
!!
!! ARGUMENTS
!!  
!! transformType - Array containing the type of transform required 
!!                 along each axis (e.g. PFFT_REAL2C or PFFT_COMPLEX)
!!
!! NOTES
!!
!! Called by gr_pfftInit.
!!
!!***

subroutine gr_pfftSpecifyTransform(transformType,baseDatType, bcTypes)
#include "constants.h"
#include "Pfft.h"
  implicit none
  integer, dimension(MDIM), intent(OUT) :: transformType
  integer, dimension(0:MDIM), intent(OUT), OPTIONAL :: baseDatType
  integer, dimension(2*MDIM), intent(IN),  OPTIONAL :: bcTypes

  !We must provide our own implementation of this function because the
  !default implementation will set 
  !transformType(IAXIS) = PFFT_REAL2C
  !transformType(JAXIS:KAXIS) = PFFT_COMPLEX
  !which will alter the lengths of certain shape arrays.
  !We use PFFT_SIN because it leaves the lengths unchanged.
  ! DEV: This should not be needed any more, since we have baseDatType!
  transformType(IAXIS:KAXIS) = PFFT_SIN

  if (present(baseDatType)) then
     baseDatType(0:MDIM) = PFFT_PCLDATA_REAL
  end if
end subroutine gr_pfftSpecifyTransform
