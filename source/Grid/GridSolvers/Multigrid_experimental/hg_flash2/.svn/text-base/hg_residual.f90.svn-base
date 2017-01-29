!******************************************************************************

! Routine:      hg_residual

! Description:  Compute residual for an approximate solution to the Poisson
!               equation.  Boundary _values_ are used in computing residuals
!               near boundaries; these have been used to arrive at the
!               solution to begin with and are either externally imposed or
!               interpolated from parent blocks (hence, no same-refinement-
!               level transfers are done).

! Parameters:   level        Level on which to compute residual
!               isrc         Source variable index
!               isoln        Solution variable index
!               ires         Index of variable to receive residual


subroutine hg_residual(level, isrc, isoln, ires)

!==============================================================================

use dBase, ONLY:  dBaseRefinementLevel, dBaseGetDataPtrSingleBlock, &
                  dBaseReleaseDataPtrSingleBlock, GC, nxb, nyb, nzb, &
                  dBaseBlockSize, ndim, nguard, dBasePropertyInteger,&
                  dBaseNodeType, maxblocks, nchild, k2d, k3d, &
                  dBaseNeighborBlockList
use dBaseDeclarations, ONLY: work

use mg_common
use perfmon
use runtime_parameters

use mpi

implicit none

integer, intent(in) :: level, isrc, isoln, ires

integer             :: b, i, j, k, lnblocks, ii, jj, kk, ierr
integer             :: nbr_blks(6)
real, pointer       :: solnData(:,:,:,:)

real                :: avg, sum, lsum, vol, lvol, nbinv, bvol, cvol, bsum

integer, parameter  :: ns = 1
integer, save       :: lrefine_min
real, save          :: cx(-ns:ns,nxb), cy(-ns:ns,nyb), cz(-ns:ns,nzb)
logical, save       :: first_call = .true.
real                :: dxinv2, dyinv2, dzinv2, bsize(3)

!==============================================================================

call timer_start("hg_residual")

if (first_call) then
  do i = 1, nxb
    cx(:,i) = (/ 1., -2., 1. /)
!    cx(:,i) = (/ -1./12., 4./3., -5./2., 4./3., -1./12. /)  ! use ns = 2
  enddo
  do j = 1, nyb
    cy(:,j) = (/ 1., -2., 1. /)
!    cy(:,j) = (/ -1./12., 4./3., -5./2., 4./3., -1./12. /)  ! use ns = 2
  enddo
  do k = 1, nzb
    cz(:,k) = (/ 1., -2., 1. /)
!    cz(:,k) = (/ -1./12., 4./3., -5./2., 4./3., -1./12. /)  ! use ns = 2
  enddo
  call get_parm_from_context("lrefine_min", lrefine_min)
  first_call = .false.
endif

lnblocks = dBasePropertyInteger("LocalNumberOfBlocks")

! Use (EXCHANGE_WORK, CONTINUE_SERIES) under the assumption that we're
! calling hg_residual immediately after a call to hg_solve_level --
! performance savings by not filling work() again

call mg_bndry(level, isoln, nguard, 0, EXCHANGE_WORK, CONTINUE_SERIES, 0)

do b = 1, lnblocks
  if (dBaseRefinementLevel(b) == level) then

    solnData => dBaseGetDataPtrSingleBlock(b, GC)

    bsize = dBaseBlockSize(b)

    do k = nguard*k3d+1, nguard*k3d+nzb
      do j = nguard*k2d+1, nguard*k2d+nyb
        do i = nguard+1, nguard+nxb
          solnData(ires,i,j,k) = solnData(isrc,i,j,k)
        enddo
      enddo
    enddo

    if (ndim == 1) then

      dxinv2 = (nxb / bsize(1))**2

      do i = nguard+1, nguard+nxb
        do ii = -ns, ns
          solnData(ires,i,1,1) = &
             solnData(ires,i,1,1) - dxinv2*cx(ii,i-nguard)*work(i+ii,1,1,b)
        enddo
      enddo

    else if (ndim == 2) then

      dxinv2 = (nxb / bsize(1))**2
      dyinv2 = (nyb / bsize(2))**2

      do j = nguard+1, nguard+nyb
        do i = nguard+1, nguard+nxb
          do ii = -ns, ns
            solnData(ires,i,j,1) = &
               solnData(ires,i,j,1) - dxinv2*cx(ii,i-nguard)*work(i+ii,j,1,b)
          enddo
          do jj = -ns, ns
            solnData(ires,i,j,1) = &
               solnData(ires,i,j,1) - dyinv2*cy(jj,j-nguard)*work(i,j+jj,1,b)
          enddo
        enddo
      enddo

    else ! ndim == 3

      dxinv2 = (nxb / bsize(1))**2
      dyinv2 = (nyb / bsize(2))**2
      dzinv2 = (nzb / bsize(3))**2

      do k = nguard+1, nguard+nzb
        do j = nguard+1, nguard+nyb
          do i = nguard+1, nguard+nxb
            do ii = -ns, ns
              solnData(ires,i,j,k) = &
                solnData(ires,i,j,k) - dxinv2*cx(ii,i-nguard)*work(i+ii,j,k,b)
            enddo
            do jj = -ns, ns
              solnData(ires,i,j,k) = &
                solnData(ires,i,j,k) - dyinv2*cy(jj,j-nguard)*work(i,j+jj,k,b)
            enddo
            do kk = -ns, ns
              solnData(ires,i,j,k) = &
                solnData(ires,i,j,k) - dzinv2*cz(kk,k-nguard)*work(i,j,k+kk,b)
            enddo
          enddo
        enddo
      enddo

    endif

    call dBaseReleaseDataPtrSingleBlock(b, solnData)

  endif
enddo

if ( (mg_bnd_cond == MG_BND_PERIODIC) ) then

  lvol = 0.
  lsum = 0.
  nbinv = 1. / real(nxb)
  if (ndim >= 2) nbinv = nbinv / real(nyb)
  if (ndim == 3) nbinv = nbinv / real(nzb)

  do b = 1, lnblocks
    if ((dBaseRefinementLevel(b) == level) .or. &
        ((dBaseNodeType(b) == 1) .and. (dBaseRefinementLevel(b) < level))) then
      bsize = dBaseBlockSize(b)
      bvol = bsize(1)
      if (ndim >= 2) bvol = bvol * bsize(2)
      if (ndim == 3) bvol = bvol * bsize(3)
      cvol = bvol * nbinv
      lvol = lvol + bvol
      solnData => dBaseGetDataPtrSingleBlock(b, GC)
      bsum = 0.
      do k = kli, kui
        do j = jli, jui
          do i = ili, iui
            bsum = bsum + solnData(ires,i,j,k)
          enddo
        enddo
      enddo
      call dBaseReleaseDataPtrSingleBlock(b, solnData)
      lsum = lsum + bsum * cvol
    endif
  enddo

  call mpi_allreduce ( lsum, sum, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr )
  call mpi_allreduce ( lvol, vol, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr )

  avg = sum / vol

  do b = 1, lnblocks
    if (dBaseRefinementLevel(b) == level) then
      solnData => dBaseGetDataPtrSingleBlock(b, GC)
      solnData(ires,:,:,:) = solnData(ires,:,:,:) - avg
      call dBaseReleaseDataPtrSingleBlock(b, solnData)
    endif
  enddo

endif

call timer_stop("hg_residual")

!==============================================================================

return
end subroutine hg_residual
