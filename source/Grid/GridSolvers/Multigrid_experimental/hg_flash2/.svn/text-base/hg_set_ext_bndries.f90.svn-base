!******************************************************************************

! Routine:      hg_set_ext_bndries

! Description:  Set external boundary values on the coarsest grid level for
!               the Huang & Greengard Poisson solver.

!               Currently only knows about given-value/Dirichlet-type boundary
!               conditions.  Also, only makes one pass, so won't work with
!               multiple coarsest-level blocks.

! Parameters:   isoln       Variable index for solution
!               level       Level to set boundary values on, or zero for all
!                           levels
!               bndtype     If 0, do regardless of boundary condition;
!                           if 1, do only if Dirichlet/given value boundaries


subroutine hg_set_ext_bndries(level, isoln, bndtype)

!==============================================================================

use dBase, ONLY:  dBaseRefinementLevel, dBaseGetDataPtrSingleBlock, &
                  dBaseReleaseDataPtrSingleBlock, GC, nxb, nyb, nzb, &
                  ndim, nguard, k2d, k3d, dBaseNeighborBlockList, &
                  dBasePropertyInteger

use dBaseIncludes, ONLY: bnd_box
use mg_common

implicit none

integer, intent(in) :: isoln, level, bndtype

integer             :: b, i, j, k, lnblocks
integer             :: nbr_blks(6)
real, pointer       :: solnData(:,:,:,:)

!==============================================================================

lnblocks = dBasePropertyInteger("LocalNumberOfBlocks")

if ((bndtype == 1) .and. (mg_bnd_cond /= MG_BND_DIRICHLET) .and. (mg_bnd_cond /= MG_BND_GIVENVAL)) return

do b = 1, lnblocks
  if ((level == 0) .or. (dBaseRefinementLevel(b) == level)) then

    nbr_blks = dBaseNeighborBlockList(b)

    solnData => dBaseGetDataPtrSingleBlock(b, GC)

    if (nbr_blks(1) <= -20) solnData(isoln,nguard,:,:) = 0.
    if (nbr_blks(2) <= -20) solnData(isoln,nguard+nxb+1,:,:) = 0.

    if (ndim >= 2) then
      if (nbr_blks(3) <= -20) solnData(isoln,:,nguard,:) = 0.
      if (nbr_blks(4) <= -20) solnData(isoln,:,nguard+nyb+1,:) = 0.
    endif

    if (ndim == 3) then
      if (nbr_blks(5) <= -20) solnData(isoln,:,:,nguard) = 0.
      if (nbr_blks(6) <= -20) solnData(isoln,:,:,nguard+nzb+1) = 0.
    endif

    call dBaseReleaseDataPtrSingleBlock(b, solnData)

  endif
enddo

!==============================================================================

return
end subroutine hg_set_ext_bndries
