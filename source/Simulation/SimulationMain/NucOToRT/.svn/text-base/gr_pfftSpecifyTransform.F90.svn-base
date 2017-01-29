!!****if* source/Simulation/SimulationMain/NucOToRT/gr_pfftSpecifyTransform
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
#include "Flash.h"
#include "constants.h"
#include "Pfft.h"
  implicit none
  integer, dimension(MDIM), intent(OUT) :: transformType
  integer, dimension(0:MDIM), intent(OUT), OPTIONAL :: baseDatType
  integer, dimension(2*MDIM), intent(IN),  OPTIONAL :: bcTypes

  !! This is temporary, assumes we are only working with 
  !! periodic boundary conditions
  if (NDIM==1) then
     transformType(IAXIS) = PFFT_REAL2C_STUFF
  else
     transformType(IAXIS) = PFFT_REAL2C
  end if
  transformType(JAXIS:KAXIS) = PFFT_COMPLEX

  if (present(baseDatType)) then
     baseDatType(0) = PFFT_PCLDATA_REAL
     baseDatType(1:MDIM) = PFFT_PCLDATA_COMPLEX
  end if
end subroutine gr_pfftSpecifyTransform
