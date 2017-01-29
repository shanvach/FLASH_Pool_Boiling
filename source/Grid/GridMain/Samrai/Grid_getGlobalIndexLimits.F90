!!****if* source/Grid/GridMain/Samrai/Grid_getGlobalIndexLimits
!!
!! NAME
!!  Grid_getGlobalIndexLimits
!!
!! SYNOPSIS
!!
!!  Grid_getGlobalIndexLimits(integer :: globalIndexLimits(MDIM))
!!  
!! DESCRIPTION 
!!
!!  Get dimensions of grid across all processors
!!
!!***


subroutine Grid_getGlobalIndexLimits(globalIndexLimits)

  use Grid_data, ONLY : gr_indexDomainMinMax

  implicit none

#include "constants.h"

  integer, dimension(MDIM), intent(OUT) :: globalIndexLimits

  globalIndexLimits(IAXIS) = gr_indexDomainMinMax(2)
  globalIndexLimits(JAXIS) = gr_indexDomainMinMax(4)
  globalIndexLimits(KAXIS) = gr_indexDomainMinMax(6)


end subroutine Grid_getGlobalIndexLimits
