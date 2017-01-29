!*******************************************************************************

!  Routine:     mg_norm()

!  Description: Compute the L2 norm of a multigrid variable on a particular
!               level.  Assume data on this level is good (e.g., because of
!               a call to mg_restrict()).  If level == 0, compute norm on
!               all levels.  If leaf_only /= 0, only compute norm on leaf
!               nodes.  If ivar == -1, compute norm on the work array.


subroutine mg_norm (level, ivar, norm, leaf_only)

!===============================================================================

use mg_common, ONLY: nodetype_save, ili, jli, kli, iui, jui, kui
use dBase, ONLY:               &
     &     ndim, nxb, nyb, nzb,      & 
     &     dBaseGetDataPtrAllBlocks, &
     &     dBaseNodeType,            &
     &     dBaseRefinementLevel,     &
     &     dBaseBlockSize,           &
     &     dBasePropertyInteger

!temporary until we can think of better way
use dBaseDeclarations, ONLY: work

use perfmon

implicit none
include "mpif.h"

integer :: level, ivar, leaf_only
real    :: norm

integer :: lb, i, j, k, ierr, lnblocks
real    :: lvol, vol, lsum, bsum, sum
real    :: nbinv, cvol, bvol
logical :: include_in_sum

integer, parameter :: MAXDIMS = 3
real, dimension(MAXDIMS) :: size

real, pointer, dimension(:,:,:,:,:), save :: unk
logical, save :: first_call = .true.

!===============================================================================

call timer_start("mg_norm")

lvol = 0.
lsum = 0.
nbinv = 1. / (real(nxb)*real(nyb)*real(nzb))

if (first_call) then
   unk => dBaseGetDataPtrAllBlocks()
   first_call = .false.
end if 

lnblocks = dBasePropertyInteger("LocalNumberOfBlocks")

do lb = 1, lnblocks
  include_in_sum = (dBaseRefinementLevel(lb) == level) .or. (level == 0)
  if (leaf_only /= 0) & 
    include_in_sum = include_in_sum .and. (nodetype_save(lb) == 1)
  if (include_in_sum) then
    size = dBaseBlockSize(lb)
    bvol = size(1)
    if (ndim >= 2) bvol = bvol * size(2)
    if (ndim == 3) bvol = bvol * size(3)
    cvol = bvol * nbinv
    lvol = lvol + bvol
    bsum = 0.
    if (ivar >= 0) then
      do k = kli, kui
        do j = jli, jui
          do i = ili, iui
            bsum = bsum + unk(ivar,i,j,k,lb)**2
          enddo
        enddo
      enddo
    else
      do k = kli, kui
        do j = jli, jui
          do i = ili, iui
            bsum = bsum + work(i,j,k,lb)**2
          enddo
        enddo
      enddo
    endif
    lsum = lsum + bsum * cvol
  endif
enddo

call mpi_allreduce ( lsum, sum, 1, MPI_DOUBLE_PRECISION, & 
                     MPI_SUM, MPI_COMM_WORLD, ierr )
!call mpi_allreduce ( lvol, vol, 1, MPI_DOUBLE_PRECISION,
!                     MPI_SUM, MPI_COMM_WORLD, ierr )

norm = sqrt(sum)  ! definition in Briggs et al., _A Multigrid Tutorial_

call timer_stop("mg_norm")

!===============================================================================

return
end
