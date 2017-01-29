!!***f* source/Grid/common/paramesh/Grid_markRefineVarThreshold
!!
!! NAME
!!  Grid_markRefineVarThreshold
!!
!! SYNOPSIS
!!  Grid_markRefineVarThreshold(real : VarVect(:, :, :, MAXBLOCKS),
!!                         real : var_th, 
!!                         integer : icmp, 
!!                         integer : lref )
!!
!! PURPOSE
!!  Refine all blocks for which a given variable (VarVect) exceeds or falls
!!  below some threshold (var_th).  The direction of the threshold is
!!  controlled by the parameter icmp.  Either blocks are brought
!!  up to a specific level of refinement or each block is refined once.
!!
!! ARGUMENTS
!!  VarVect(nxb, nyb, nzb, maxblocks)   the variable of interest
!!
!!  var_th                              the limit on the variable
!! 
!!  icmp   icmp < 0  refine if the variable is less than var_th
!!         icmp >= 0 refine if the variable is greater then var_th
!! 
!!  lref   max refinement level
!!

subroutine Grid_markRefineVarThreshold (VarVect, var_th, icmp,lref)

  implicit none
#include "constants.h"
! Arguments

  real, dimension(NXB, NYB, NZB, MAXBLOCKS), intent(IN) :: VarVect
  real,    intent(IN) :: var_th
  integer, intent(IN) :: icmp, lref

  return
end subroutine Grid_markRefineVarThreshold
