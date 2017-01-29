!*******************************************************************************

! Routine:      mg_guardcell

! Description:  Version of amr_guardcell() for the multigrid solver.
!               The only differences are (1) it only operates on the
!               work array and (2) it calls the special routine
!               mg_set_ext_bndry() to handle exterior guard cells.

!               See documentation for amr_guardcell() for more information.


subroutine mg_guardcell(mype, ivar, nlayers, simtime, idiag, idir, extrap)

!===============================================================================

use dBase, ONLY: dBaseGetDataPtrAllBlocks, dBasePropertyInteger,     &
                 dBasePropertyReal, dBaseRefinementLevel,            &
                 dBaseTreePtrNodeType, dBaseNeighborBlockList,       &
                 nguard_work, maxblocks_tr, ndim, k2d, k3d, nxb, nyb, nguard

use mg_common
use perfmon
use dBaseDeclarations, only: work
use runtime_scratch
implicit none

include 'mpif.h'

integer :: indexPlotNumber, PTNumber, ierr
integer :: idebug
integer :: mype, nlayers, idiag, idir, extrap
integer ivar
integer ::  lnblocks, lb, i,j,k
real    :: simtime
real,    pointer, dimension(:,:,:,:,:), save :: solnData
logical, save :: first_call = .true.
!===============================================================================

if (first_call) then
   solnData => dBaseGetDataPtrAllBlocks()
   first_call = .false.
end if

!call timer_start("in_mg_guardcell")
if (nlayers > nguard_work) &
  call abort_flash("mg_guardcell: you must set nlayers <= nguard_work")

! Put recognizably bad data in corners to permit recognition of
! corner guard cells diagonally opposite coarser blocks.

if (idiag == 1) call amr_mark_edges(mype, 2)

! Restrict data from leaf blocks to their parents. This will enable the
! parent blocks to provide guard cell data to their neighbors in cases
! where the leaf blocks have coarser neighbors.

call amr_restrict(mype, 2, 0)


! Blocks provide guard cell data to all their neighbors which share their
! level of refinement.

!call timer_start("guardcell srl")
call amr_guardcell_srl(mype, 2, nlayers, idiag, idir)
!call timer_stop("guardcell srl")

! Apply boundary conditions for c_to_f by filling guard cells at external
! boundaries with extrapolation turned on.
call mg_set_ext_bndry(idiag, idir, 1)

! Guard cell data is sent to any neighbor blocks with finer resolution.
!call timer_start("guardcell c_to_f")
call amr_guardcell_c_to_f(mype, 2, nlayers, idiag, idir)
!call timer_stop("guardcell c_to_f")
!call timer_stop("in_mg_guardcell")

! Apply boundary conditions by filling guard cells at external boundaries.
call mg_set_ext_bndry(idiag, idir, extrap)

!===============================================================================

return
end

