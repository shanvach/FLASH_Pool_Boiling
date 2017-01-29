!*******************************************************************************

! Routine:      poisson_image_mass

! Description:  Given a zero-boundary potential in the solution variable with
!               index isoln, compute the surface density of the image mass
!               required to produce this potential.  The result is written into
!               the solution variable with index iiden.


subroutine poisson_image_mass (isoln, iiden)

!===============================================================================

use mg_common
use dBase, ONLY:  dBasePropertyInteger, dBaseNeighborBlockList, &
                  dBaseNodeType, dBaseGetDataPtrSingleBlock, &
                  dBaseBlockSize, nxb, nyb, nzb, nguard, GC, &
                  dBaseReleaseDataPtrSingleBlock
implicit none

integer       :: isoln, iiden

integer       :: i, j, k, lb, lnblocks, nbrs(6), ndim
real, pointer :: solnData(:,:,:,:)
real          :: bsize(3), dx, dy, dz

!===============================================================================

lnblocks = dBasePropertyInteger("LocalNumberOfBlocks")
ndim = dBasePropertyInteger("Dimensionality")

do lb = 1, lnblocks
  solnData => dBaseGetDataPtrSingleBlock(lb, GC)
  solnData(iiden,:,:,:) = 0.
  if (dBaseNodeType(lb) == 1) then
    nbrs = dBaseNeighborBlockList(lb)
    bsize= dBaseBlockSize(lb)

    dx = bsize(1) / nxb

    if (nbrs(1) <= -20) then       ! -x boundary
      solnData(iiden,nguard,:,:) = &
                         3.5/dx * solnData(isoln,nguard+1,:,:)  &
                         -0.5/dx * solnData(isoln,nguard+2,:,:)
      solnData(iiden,nguard+1,jli:jui,kli:kui) = &
                         0.5/dx * solnData(iiden,nguard,jli:jui,kli:kui)
    endif

    if (nbrs(2) <= -20) then       ! +x boundary
      solnData(iiden,nguard+nxb+1,:,:) = &
                         3.5/dx * solnData(isoln,nguard+nxb,:,:)  &
                         -0.5/dx * solnData(isoln,nguard+nxb-1,:,:)
      solnData(iiden,nguard+nxb,jli:jui,kli:kui) = &
                         0.5/dx * solnData(iiden,nguard+nxb+1,jli:jui,kli:kui)
    endif

    if (ndim >= 2) then
      dy = bsize(2) / nyb
      if (nbrs(3) <= -20) then       ! -y boundary
        solnData(iiden,:,nguard,:) = &
                         3.5/dy * solnData(isoln,:,nguard+1,:)  &
                         -0.5/dy * solnData(isoln,:,nguard+2,:)
        solnData(iiden,ili:iui,nguard+1,kli:kui) = &
                         solnData(iiden,ili:iui,nguard+1,kli:kui) + &
                         0.5/dy * solnData(iiden,ili:iui,nguard,kli:kui)
      endif

      if (nbrs(4) <= -20) then       ! +y boundary
        solnData(iiden,:,nguard+nyb+1,:) = &
                         3.5/dy * solnData(isoln,:,nguard+nyb,:)  &
                         -0.5/dy * solnData(isoln,:,nguard+nyb-1,:)
        solnData(iiden,ili:iui,nguard+nyb,kli:kui) = &
                         solnData(iiden,ili:iui,nguard+nyb,kli:kui) + &
                         0.5/dy * solnData(iiden,ili:iui,nguard+nyb+1,kli:kui)
      endif
    endif

    if (ndim == 3) then
      dz = bsize(3) / nzb
      if (nbrs(5) <= -20) then       ! -z boundary
        solnData(iiden,:,:,nguard) = &
                         3.5/dz * solnData(isoln,:,:,nguard+1)  &
                         -0.5/dz * solnData(isoln,:,:,nguard+2)
        solnData(iiden,ili:iui,jli:jui,nguard+1) = &
                         solnData(iiden,ili:iui,jli:jui,nguard+1) + &
                         0.5/dz * solnData(iiden,ili:iui,jli:jui,nguard)
      endif

      if (nbrs(6) <= -20) then       ! +z boundary
        solnData(iiden,:,:,nguard+nzb+1) = &
                         3.5/dz * solnData(isoln,:,:,nguard+nzb)  &
                         -0.5/dz * solnData(isoln,:,:,nguard+nzb-1)
        solnData(iiden,ili:iui,jli:jui,nguard+nzb) = &
                         solnData(iiden,ili:iui,jli:jui,nguard+nzb) + &
                         0.5/dz * solnData(iiden,ili:iui,jli:jui,nguard+nzb+1)
      endif
    endif

  endif

  call dBaseReleaseDataPtrSingleBlock(lb, solnData)

enddo

!===============================================================================

return
end
