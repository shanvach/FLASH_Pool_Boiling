!******************************************************************************

! Routine:      hg_solve_level

! Description:  Invoke a specified single-block Poisson solver on all blocks
!               at a given level of refinement.  Boundary conditions for the
!               blocks are assumed to have been set already and are stored in
!               the first layer of guard sells of the solution variable.

!               Currently only knows about given-value/Dirichlet-type boundary
!               conditions.  Also, only makes one pass, so won't work with
!               multiple coarsest-level blocks.

! Parameters:   level       Level on which to solve
!               isrc        Variable index for source function
!               isoln       Variable index for solution
!               SolveBlock  A routine to call a single-block direct Poisson
!                           solver once boundary values have been set.
!               LeafFlag    Flag controlling whether leaf blocks are solved.
!                           If 0, solve all blocks on this level.  If 1, solve
!                           only leaf blocks.  If 2, solve only parent blocks.


subroutine hg_solve_level(level, isrc, isoln, SolveBlock, LeafFlag)

!==============================================================================

use mg_common

use dBase, ONLY:  dBaseRefinementLevel, dBaseGetDataPtrSingleBlock, &
                  dBaseReleaseDataPtrSingleBlock, GC, nxb, nyb, nzb, &
                  dBaseBlockSize, ndim, nguard, k2d, k3d, &
                  dBaseNeighborBlockList, dBasePropertyInteger, dBaseNodeType
use dBaseDeclarations, ONLY: work
use perfmon
use runtime_parameters

use mpi

implicit none

integer, intent(in) :: isrc, isoln, level, LeafFlag
external            :: SolveBlock

integer                      :: b, i, j, k, lnblocks, n, nbr_blks(6)
real, pointer                :: solnData(:,:,:,:)
real, dimension(nxb,nyb,nzb) :: soln
real, dimension(3)           :: block_size, zone_size
real                         :: dx, dy, dz, c, cx, cy, cz
real                         :: avg, sum, lsum, vol, lvol, nbinv, bvol, cvol, bsum
logical                      :: SolveThisBlock
integer                      :: bnd_type, ierr

logical, save                :: first_call = .true.
integer, save                :: nxl1, nxl2, nyl1, nyl2, nzl1, nzl2
integer, save                :: nxr1, nxr2, nyr1, nyr2, nzr1, nzr2
integer, save                :: lrefine_min, MyPE

!==============================================================================
call timer_start("hg_solve_level")
call timer_start("fft")

if (first_call) then
  nxl1 = 1
  nxl2 = max(nxb/2, 1)
  nxr1 = max(nxb/2+1, 1)
  nxr2 = nxb
  nyl1 = 1
  nyl2 = max(nyb/2, 1)
  nyr1 = max(nyb/2+1, 1)
  nyr2 = nyb
  nzl1 = 1
  nzl2 = max(nzb/2, 1)
  nzr1 = max(nzb/2+1, 1)
  nzr2 = nzb
  call get_parm_from_context("lrefine_min", lrefine_min)
  MyPE = dBasePropertyInteger("MyProcessor")
  first_call = .false.
endif

lnblocks = dBasePropertyInteger("LocalNumberOfBlocks")

do b = 1, lnblocks

  SolveThisBlock = (dBaseRefinementLevel(b) == level)
  if (LeafFlag == 1) then
    SolveThisBlock = (SolveThisBlock .and. (dBaseNodeType(b) == 1))
  else if (LeafFlag == 2) then
    SolveThisBlock = (SolveThisBlock .and. (dBaseNodeType(b) /= 1))
  endif

  if (SolveThisBlock) then

    block_size = dBaseBlockSize(b)
    dx = block_size(1) / nxb
    if (ndim >= 2) then
      dy = block_size(2) / nyb
    else
      dy = 1.
    endif
    if (ndim == 3) then
      dz = block_size(3) / nzb
    else
      dz = 1.
    endif

    solnData => dBaseGetDataPtrSingleBlock(b, GC)

    do k = 1, nzb
      do j = 1, nyb
        do i = 1, nxb
          soln(i,j,k) = solnData(isrc,i+nguard,j+k2d*nguard,k+k3d*nguard)
        enddo
      enddo
    enddo

    if (mg_bnd_cond == MG_BND_PERIODIC) then
      ! Only the coarsest mesh level is treated as periodic; the rest obtain
      ! boundary values by interpolation and are treated as Dirichlet.
      if (level == 1) then
        bnd_type = 0
      else
        bnd_type = 1
      endif
    else if ((mg_bnd_cond == MG_BND_DIRICHLET) .or. &
             (mg_bnd_cond == MG_BND_GIVENVAL)) then
      bnd_type = 1
    endif

    ! Adjust the source function to account for periodic boundary conditions.

!    if (bnd_type == 0) then
!
!      avg = sum(soln) / (nxb*nyb*nzb)
!      soln = soln - avg
!
!    endif

    ! Adjust the source function to account for given-value boundary
    ! conditions.

    if (bnd_type == 1) then

      nbr_blks = dBaseNeighborBlockList(b)

        soln(1,:,:)   = soln(1,:,:)   - 2.*solnData(isoln,nguard, &
                                                1+nguard*k2d:nyb+nguard*k2d, &
                                                1+nguard*k3d:nzb+nguard*k3d)/dx**2
        soln(nxb,:,:) = soln(nxb,:,:) - 2.*solnData(isoln,nguard+nxb+1, &
                                                1+nguard*k2d:nyb+nguard*k2d, &
                                                1+nguard*k3d:nzb+nguard*k3d)/dx**2

      if (ndim >= 2) then
          soln(:,1,:)   = soln(:,1,:)   - 2.*solnData(isoln,nguard+1:nguard+nxb, &
                                                nguard, &
                                                1+nguard*k3d:nzb+nguard*k3d)/dy**2
          soln(:,nyb,:) = soln(:,nyb,:) - 2.*solnData(isoln,nguard+1:nguard+nxb, &
                                                nguard+nyb+1, &
                                                1+nguard*k3d:nzb+nguard*k3d)/dy**2
      endif

      if (ndim == 3) then
          soln(:,:,1)   = soln(:,:,1)   - 2.*solnData(isoln,nguard+1:nguard+nxb, &
                                                  1+nguard:nyb+nguard, &
                                                  nguard           )/dz**2
          soln(:,:,nzb) = soln(:,:,nzb) - 2.*solnData(isoln,nguard+1:nguard+nxb, &
                                                  1+nguard:nyb+nguard, &
                                                  nguard+nzb+1           )/dz**2
      endif

    endif

    call SolveBlock (soln, nxb, nyb, nzb, dx, dy, dz, bnd_type, level)

    do k = 1, nzb
      do j = 1, nyb
        do i = 1, nxb
          solnData(isoln,i+nguard,j+k2d*nguard,k+k3d*nguard) = soln(i,j,k)
          work(i+nguard,j+k2d*nguard,k+k3d*nguard,b) = soln(i,j,k)
        enddo
      enddo
    enddo

    call dBaseReleaseDataPtrSingleBlock(b, solnData)

  endif
enddo
call timer_stop("fft")

!TEST - BOUNDARY RELAXATION
call timer_start("relaxation")

do n = 1, 2

if (n == 1) then
  call mg_bndry(level, isoln, 1, 0, EXCHANGE_WORK, BEGIN_SERIES, 0)
else
  call mg_bndry(level, isoln, 1, 0, EXCHANGE_WORK, CONTINUE_SERIES, 0)
endif

do b = 1, lnblocks

  SolveThisBlock = (dBaseRefinementLevel(b) == level)
  if (LeafFlag == 1) then
    SolveThisBlock = (SolveThisBlock .and. (dBaseNodeType(b) == 1))
  else if (LeafFlag == 2) then
    SolveThisBlock = (SolveThisBlock .and. (dBaseNodeType(b) /= 1))
  endif

  if (SolveThisBlock) then

    block_size = dBaseBlockSize(b)
    dx = block_size(1) / nxb
    if (ndim >= 2) then
      dy = block_size(2) / nyb
    else
      dy = 1.
    endif
    if (ndim == 3) then
      dz = block_size(3) / nzb
    else
      dz = 1.
    endif

    solnData => dBaseGetDataPtrSingleBlock(b, GC)

    if (ndim == 1) then

      cx = 0.5
      c  = 0.5 * dx**2

      do i = nguard+1, nguard+2
        work(i,1,1,b) = cx*(work(i-1,1,1,b) + work(i+1,1,1,b)) - &
                        c*solnData(isrc,i,1,1)
        solnData(isoln,i,1,1) = work(i,1,1,b)
      enddo

      do i = nguard+nxb-1, nguard+nxb
        work(i,1,1,b) = cx*(work(i-1,1,1,b) + work(i+1,1,1,b)) - &
                        c*solnData(isrc,i,1,1)
        solnData(isoln,i,1,1) = work(i,1,1,b)
      enddo

    else if (ndim == 2) then

      c  = 0.5 / (dx**2 + dy**2)
      cx = dy**2 * c
      cy = dx**2 * c
      c  = dx**2 * dy**2 * c

      do j = nguard+1, nguard+2
        do i = nguard+1, nguard+nxb
          work(i,j,1,b) = &
            cx*(work(i-1,j,1,b) + work(i+1,j,1,b)) + &
            cy*(work(i,j-1,1,b) + work(i,j+1,1,b)) - &
            c*solnData(isrc,i,j,1)
          solnData(isoln,i,j,1) = work(i,j,1,b)
        enddo
      enddo
      do j = nguard+nyb-1, nguard+nyb
        do i = nguard+1, nguard+nxb
          work(i,j,1,b) = &
            cx*(work(i-1,j,1,b) + work(i+1,j,1,b)) + &
            cy*(work(i,j-1,1,b) + work(i,j+1,1,b)) - &
            c*solnData(isrc,i,j,1)
          solnData(isoln,i,j,1) = work(i,j,1,b)
        enddo
      enddo
      do j = nguard+1, nguard+nyb
        do i = nguard+1, nguard+2
          work(i,j,1,b) = &
            cx*(work(i-1,j,1,b) + work(i+1,j,1,b)) + &
            cy*(work(i,j-1,1,b) + work(i,j+1,1,b)) - &
            c*solnData(isrc,i,j,1)
          solnData(isoln,i,j,1) = work(i,j,1,b)
        enddo
        do i = nguard+nxb-1, nguard+nxb
          work(i,j,1,b) = &
            cx*(work(i-1,j,1,b) + work(i+1,j,1,b)) + &
            cy*(work(i,j-1,1,b) + work(i,j+1,1,b)) - &
            c*solnData(isrc,i,j,1)
          solnData(isoln,i,j,1) = work(i,j,1,b)
        enddo
      enddo

    else ! ndim == 3

      c  = 0.5 / (dy**2*dz**2 + dx**2*dz**2 + dx**2*dy**2)
      cx = dy**2*dz**2 * c
      cy = dx**2*dz**2 * c
      cz = dx**2*dy**2 * c
      c  = dx**2*dy**2*dz**2 * c

      do k = nguard+1, nguard+2
        do j = nguard+1, nguard+nyb
          do i = nguard+1, nguard+nxb
            work(i,j,k,b) = &
              cx*(work(i-1,j,k,b) + work(i+1,j,k,b)) + &
              cy*(work(i,j-1,k,b) + work(i,j+1,k,b)) + &
              cz*(work(i,j,k-1,b) + work(i,j,k+1,b)) - &
              c*solnData(isrc,i,j,k)
            solnData(isoln,i,j,k) = work(i,j,k,b)
          enddo
        enddo
      enddo
      do k = nguard+nzb-1, nguard+nzb
        do j = nguard+1, nguard+nyb
          do i = nguard+1, nguard+nxb
            work(i,j,k,b) = &
              cx*(work(i-1,j,k,b) + work(i+1,j,k,b)) + &
              cy*(work(i,j-1,k,b) + work(i,j+1,k,b)) + &
              cz*(work(i,j,k-1,b) + work(i,j,k+1,b)) - &
              c*solnData(isrc,i,j,k)
            solnData(isoln,i,j,k) = work(i,j,k,b)
          enddo
        enddo
      enddo
      do k = nguard+1, nguard+nzb
        do j = nguard+1, nguard+2
          do i = nguard+1, nguard+nxb
            work(i,j,k,b) = &
              cx*(work(i-1,j,k,b) + work(i+1,j,k,b)) + &
              cy*(work(i,j-1,k,b) + work(i,j+1,k,b)) + &
              cz*(work(i,j,k-1,b) + work(i,j,k+1,b)) - &
              c*solnData(isrc,i,j,k)
            solnData(isoln,i,j,k) = work(i,j,k,b)
          enddo
        enddo
        do j = nguard+nyb-1, nguard+nyb
          do i = nguard+1, nguard+nxb
            work(i,j,k,b) = &
              cx*(work(i-1,j,k,b) + work(i+1,j,k,b)) + &
              cy*(work(i,j-1,k,b) + work(i,j+1,k,b)) + &
              cz*(work(i,j,k-1,b) + work(i,j,k+1,b)) - &
              c*solnData(isrc,i,j,k)
            solnData(isoln,i,j,k) = work(i,j,k,b)
          enddo
        enddo
        do j = nguard+1, nguard+nyb
          do i = nguard+1, nguard+2
            work(i,j,k,b) = &
              cx*(work(i-1,j,k,b) + work(i+1,j,k,b)) + &
              cy*(work(i,j-1,k,b) + work(i,j+1,k,b)) + &
              cz*(work(i,j,k-1,b) + work(i,j,k+1,b)) - &
              c*solnData(isrc,i,j,k)
            solnData(isoln,i,j,k) = work(i,j,k,b)
          enddo
          do i = nguard+nxb-1, nguard+nxb
            work(i,j,k,b) = &
              cx*(work(i-1,j,k,b) + work(i+1,j,k,b)) + &
              cy*(work(i,j-1,k,b) + work(i,j+1,k,b)) + &
              cz*(work(i,j,k-1,b) + work(i,j,k+1,b)) - &
              c*solnData(isrc,i,j,k)
            solnData(isoln,i,j,k) = work(i,j,k,b)
          enddo
        enddo
      enddo

    endif

    call dBaseReleaseDataPtrSingleBlock(b, solnData)

  endif
enddo
enddo
call timer_stop("relaxation")
!TEST - BOUNDARY RELAXATION
call timer_stop("hg_solve_level")

!TEST - NORMALIZATION
if ( (mg_bnd_cond == MG_BND_PERIODIC) ) then

  lvol = 0.
  lsum = 0.
  nbinv = 1. / real(nxb)
  if (ndim >= 2) nbinv = nbinv / real(nyb)
  if (ndim == 3) nbinv = nbinv / real(nzb)

  do b = 1, lnblocks
    if ((dBaseRefinementLevel(b) == level) .or. &
        ((dBaseNodeType(b) == 1) .and. (dBaseRefinementLevel(b) < level))) then
      block_size = dBaseBlockSize(b)
      bvol = block_size(1)
      if (ndim >= 2) bvol = bvol * block_size(2)
      if (ndim == 3) bvol = bvol * block_size(3)
      cvol = bvol * nbinv
      lvol = lvol + bvol
      bsum = 0.
      do k = kli, kui
        do j = jli, jui
          do i = ili, iui
            bsum = bsum + work(i,j,k,b)
          enddo
        enddo
      enddo
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
      work(:,:,:,b) = work(:,:,:,b) - avg
      solnData => dBaseGetDataPtrSingleBlock(b, GC)
      solnData(isoln,:,:,:) = solnData(isoln,:,:,:) - avg
      call dBaseReleaseDataPtrSingleBlock(b, solnData)
    endif
  enddo

endif
!TEST - NORMALIZATION

!==============================================================================

return
end subroutine hg_solve_level
