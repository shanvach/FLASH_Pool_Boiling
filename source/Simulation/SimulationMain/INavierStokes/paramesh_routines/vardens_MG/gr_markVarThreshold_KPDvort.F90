!!****if* source/Simulation/SimulationMain/INavierStokes/paramesh_routines/vardens_MG/gr_markVarThreshold_KPD
!!
!! NAME
!!  gr_markVarThreshold
!!
!! SYNOPSIS
!!  gr_markVarThreshold(  integer(in) :: Var
!!                         real(in)   :: var_th, 
!!                         integer(in) :: icmp, 
!!                         integer(in) :: lref )
!!
!! PURPOSE
!!  Refine all blocks for which a given variable (Var) exceeds or falls
!!  below some threshold (var_th).  The direction of the threshold is
!!  controlled by the parameter icmp.  Either blocks are brought
!!  up to a specific level of refinement or each block is refined once.
!!
!! ARGUMENTS
!!  Var -    the variable of interest
!!
!!  var_th  -     the limit on the variable
!! 
!!  icmp  -  icmp < 0  refine if the variable is less than var_th
!!         icmp >= 0 refine if the variable is greater then var_th
!! 
!!   lref -       If > 0, bring all qualifying blocks to this level of refinement.
!!
!!               If <= 0, refine qualifying blocks once.
!!
!! NOTES
!! 
!!  This routine has not yet been tested and should be used only as a guideline for
!!  a user's implementation.
!!
!!***

!!REORDER(5):unk

subroutine gr_markVarThreshold_KPDvort (Var, var_th, var_th2,lref)

!                                                                                          Added KPD
  use tree, ONLY : refine, derefine, lrefine, lnblocks, nodetype,lrefine_max, lrefine_min, grid_changed,coord

  use Driver_data, ONLY: dr_nstep

  use physicaldata, ONLY : unk

  implicit none

#include "constants.h"
#include "Flash.h"

! Arguments

  integer, intent(IN) :: Var
  real,    intent(IN) :: var_th, var_th2
  integer, intent(IN) :: lref

! Local data

  integer :: b,icount,zStart,zEnd
  logical :: Grid_mark
  real :: kpd_thresh

  real, dimension(MDIM) :: blockCenter

!-------------------------------------------------------------------------------

! Loop over all leaf-node blocks.

  icount = 0

  do b = 1, lnblocks
    if (nodetype(b) == LEAF) then

      icount = icount+1

      blockCenter = coord(:,b)

      !print*,"Ref Block Info: ",b,blockCenter(1),blockCenter(2),blockCenter(3)

      if (NZB .eq. 1) then
         zStart = 1
         zEnd   = 1
      else
         zStart = NGUARD+1
         zEnd   = NZB+NGUARD
      end if

if ( blockCenter(2) .LT. 0.0 .AND. &
(maxval(unk(var,NGUARD+1:NXB+NGUARD,NGUARD+1:NYB+NGUARD,NGUARD+1:NZB+NGUARD,b)) > var_th .OR. &
  minval(unk(var,NGUARD+1:NXB+NGUARD,NGUARD+1:NYB+NGUARD,NGUARD+1:NZB+NGUARD,b)) < var_th2) ) then

        if (lrefine(b) < lref ) then
          refine(b)   = .true.
          derefine(b) = .false.
        else if (lrefine(b) == lref) then
          derefine(b) = .false.
        else if (lref <= 0) then
          refine(b) = .true.
        endif

      else !kpd - Block Derefinement

        if (lrefine(b) > lrefine_min) then            !- kpd - New Derefinement
        !  refine(b)   = .false.                       !- kpd - New Derefinement
          derefine(b) = .true.                        !- kpd - New Derefinement
          !KPD - This doesn't ALWAYS mean that the grid WILL BE DEREFINED!!! Neighbor 2 rule...
        !else if (lrefine(b) <= lrefine_min) then      !- kpd - New Derefinement
        !  derefine(b) = .false.                       !- kpd - New Derefinement
        end if

      endif

    endif
  enddo

!-------------------------------------------------------------------------------

  return
end subroutine gr_markVarThreshold_KPDvort
