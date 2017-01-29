!*******************************************************************************

!  Routine:     hg_init()

!  Description: Huang-Greengard solver multigrid initialization routine.


subroutine hg_init

!===============================================================================

use mg_common
use runtime_parameters
use perfmon
use mpi

use dBase, ONLY: nxb, nyb, nzb, nguard, k2d, k3d, ndim, dBaseNodeType, &
                 dBaseRefinementLevel, dBasePropertyInteger, &
                 dBaseNeighborBlockList, dBaseGetDataPtrAllBlocks, &
                 dBaseTreePtrNodetype, dBaseTreePtrNewChild

use dBaseDeclarations, ONLY: geom_cartesian, geom_planar, geom_cylrad, &
                             geom_sphrad, geom_cylang, geom_sphtheta, &
                             geom_sphphi, lrefine

implicit none

integer                         :: lb, i, ierr, lnblocks
integer                         :: mylrefmin, mylrefmax
integer, parameter              :: MAXDIM = 3
integer                         :: nbr_blks(2*MAXDIM)

logical, save                   :: FirstCall = .true.
integer                         :: igeomx, igeomy, igeomz
integer                         :: nblockx, nblocky, nblockz
integer, pointer, dimension(:)  :: nodetype
logical, pointer, dimension(:)  :: newchild

!===============================================================================

! Save the grid structure coming into the algorithm.

lnblocks = dBasePropertyInteger("LocalNumberOfBlocks")
nodetype => dBaseTreePtrNodeType()
newchild => dBaseTreePtrNewChild()

do lb = 1, lnblocks
  nodetype_save(lb) = nodetype(lb)
  newchild_save(lb) = newchild(lb)
end do

! Determine minimum and maximum levels of refinement in the mesh.

mylrefmin = minval(lrefine(1:lnblocks))
call mpi_allreduce(mylrefmin, mesh_lrefmin, 1, MPI_INTEGER, MPI_MIN, &
                   MPI_COMM_WORLD, ierr)

mylrefmax = maxval(lrefine(1:lnblocks))
call mpi_allreduce(mylrefmax, mesh_lrefmax, 1, MPI_INTEGER, MPI_MAX, &
                   MPI_COMM_WORLD, ierr)

!report on current maximum refinement level
if (dBasePropertyInteger("MyProcessor")==0) &
   print *, 'hg_init: max refine level = ', mesh_lrefmax

! Determine index ranges for interior zones.

ili = 1 + nguard
iui = nxb + nguard
jli = 1 + k2d*nguard
jui = nyb + k2d*nguard
kli = 1 + k3d*nguard
kui = nzb + k3d*nguard

! Determine index ranges for exterior zones.

ile = 1
iue = nxb + 2*nguard
jle = 1
jue = nyb + 2*nguard*k2d
kle = 1
kue = nzb + 2*nguard*k3d

! Determine mesh geometry and decide whether we support it (only on the first
! call).

if (FirstCall) then
  call get_parm_from_context("igeomx", igeomx)
  call get_parm_from_context("igeomy", igeomy)
  call get_parm_from_context("igeomz", igeomz)
  call get_parm_from_context("quadrant", quadrant)
  select case (ndim)
  case (1)
    if (igeomx == geom_cartesian) then
      mg_geometry = MG_GEOM_1DCARTESIAN
!    elseif (igeomx == geom_cylrad) then
!      mg_geometry = MG_GEOM_1DCYLINDRICAL
!    elseif (igeomx == geom_sphrad) then
!      mg_geometry = MG_GEOM_1DSPHERICAL
    else
      mg_geometry = MG_GEOM_INVALID
    endif
  case (2)
    if ((igeomx == geom_cartesian) .and. &
        (igeomy == geom_cartesian)) then
      mg_geometry = MG_GEOM_2DCARTESIAN
!    elseif ((igeomx == geom_cylrad) .and. &
!            (igeomy == geom_planar)) then
!      mg_geometry = MG_GEOM_2DCYLAXISYM
!    elseif ((igeomx == geom_cylrad) .and. &
!            (igeomy == geom_cylang)) then
!      mg_geometry = MG_GEOM_2DCYLPOLAR
!    elseif ((igeomx == geom_sphrad) .and. &
!            (igeomy == geom_sphtheta)) then
!      mg_geometry = MG_GEOM_2DSPHAXISYM
    else
      mg_geometry = MG_GEOM_INVALID
    endif
  case (3)
    if ((igeomx == geom_cartesian) .and. &
        (igeomy == geom_cartesian) .and. &
        (igeomz == geom_cartesian)) then
      mg_geometry = MG_GEOM_3DCARTESIAN
!    elseif ((igeomx == geom_cylrad) .and. &
!            (igeomy == geom_cylang) .and. &
!            (igeomz == geom_planar)) then
!      mg_geometry = MG_GEOM_3DCYLINDRICAL
!    elseif ((igeomx == geom_sphrad) .and. &
!            (igeomy == geom_sphtheta) .and. &
!            (igeomz == geom_sphphi)) then
!      mg_geometry = MG_GEOM_3DSPHERICAL
    else
      mg_geometry = MG_GEOM_INVALID
    endif
  end select
  if (mg_geometry == MG_GEOM_INVALID) then
    call abort_flash ('hg_init:  unsupported grid geometry')
  endif

! Make sure we only have one mesh block on the coarsest level.

  call get_parm_from_context("Nblockx", nblockx)
  call get_parm_from_context("Nblocky", nblocky)
  call get_parm_from_context("Nblockz", nblockz)

  if ((nblockx /= 1) .or. (nblocky /= 1) .or. (nblockz /= 1)) then
    call abort_flash("hg_init:  only one block allowed on coarsest level")
  endif

  FirstCall = .false.
endif

!===============================================================================

return
end subroutine hg_init
