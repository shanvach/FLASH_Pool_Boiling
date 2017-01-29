!*******************************************************************************

!  Routine:     hg_init_src()

!  Description: Initializes the right-hand side of the equation to be solved
!               with the Huang & Greengard multigrid solver.
!
!               Note that the source scaling and mean subtraction operations
!               are performed directly on the source variable, so when using
!               periodic or given-value boundaries you should pass the source
!               into hg_solve() in a temporary mesh variable.


subroutine hg_init_src (isrc, isoln)

!===============================================================================

use mg_common
use physical_constants
use runtime_parameters
use dBase, ONLY:                        &
     &     nxb, nyb, nzb,  ndim, nguard,&
     &     dBasePropertyInteger,        &
     &     dBaseGetDataPtrAllBlocks,    &
     &     dBaseNodeType, dBaseBlockSize, &
           dBaseNeighborBlockList,      &
           dBaseRefinementLevel
use mpi

use runtime_scratch

implicit none

integer, intent(in) :: isrc, isoln

integer       :: lb, i, j, k, ierr, lnblocks
real          :: lvol, vol, lsum, bsum, sum
real          :: nbinv, cvol, bvol, delx, dely, delz

integer, parameter :: MAXDIMS = 3
real               :: size(MAXDIMS), nbrs(2*MAXDIMS)

real, pointer, dimension(:,:,:,:,:) :: solnData

!===============================================================================

solnData => dBaseGetDataPtrAllBlocks()

lnblocks = dBasePropertyInteger("LocalNumberOfBlocks")

! Need to restrict the source so it is valid at all levels.

do i = mesh_lrefmax, 2, -1
  call hg_restrict(i, isrc, isrc)
enddo

! For periodic or Neumann boundary conditions, compute the average value of
! the source and subtract it from the source.

src_avg = 0.

if ((mg_bnd_cond == MG_BND_PERIODIC) .or. &
    (mg_bnd_cond == MG_BND_NEUMANN)) then

  lvol = 0.
  lsum = 0.
  nbinv = 1. / real(nxb)
  if (ndim >= 2) nbinv = nbinv / real(nyb)
  if (ndim == 3) nbinv = nbinv / real(nzb)

  do lb = 1, lnblocks
!    if (dBaseNodeType(lb) == 1) then
    if (dBaseRefinementLevel(lb) == 1) then
      size = dBaseBlockSize(lb)
      bvol = size(1)
      if (ndim >= 2) bvol = bvol * size(2)
      if (ndim == 3) bvol = bvol * size(3)
      cvol = bvol * nbinv
      lvol = lvol + bvol
      bsum = 0.
      do k = kli, kui
        do j = jli, jui
          do i = ili, iui
            bsum = bsum + solnData(isrc,i,j,k,lb)
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

  src_avg = sum / vol

  do lb = 1, lnblocks
    solnData(isrc,:,:,:,lb) = solnData(isrc,:,:,:,lb) - src_avg
  enddo

! For given-value boundary conditions, we must subtract exterior boundary
! values (assumed to be stored in the first layer of boundary zones for
! the solution variable).

else if (mg_bnd_cond == MG_BND_GIVENVAL) then

  do lb = 1, lnblocks

    nbrs = dBaseNeighborBlockList(lb)
    size = dBaseBlockSize(lb)
    delx = size(1) / nxb
    if (ndim >= 2) dely = size(2) / nyb
    if (ndim == 3) delz = size(3) / nzb

    if ((mg_geometry == MG_GEOM_1DCARTESIAN) .or. &
        (mg_geometry == MG_GEOM_2DCARTESIAN) .or. &
        (mg_geometry == MG_GEOM_3DCARTESIAN)) then

      if (nbrs(1) <= -20) then
        solnData(isrc,nguard+1,:,:,lb) = &
          solnData(isrc,nguard+1,:,:,lb) - &
          2.*solnData(isoln,nguard,:,:,lb)/delx**2
      endif

    endif

    if ((mg_geometry /= MG_GEOM_2DCYLAXISYM) .or. &
        (.not. quadrant)) then

      if ((ndim >= 2) .and. (nbrs(3) <= -20)) then
        solnData(isrc,:,nguard+1,:,lb) = &
          solnData(isrc,:,nguard+1,:,lb) - &
          2.*solnData(isoln,:,nguard,:,lb)/dely**2
      endif

    endif

    if ((ndim == 3) .and. (nbrs(5) <= -20)) then
      solnData(isrc,:,:,nguard+1,lb) = &
        solnData(isrc,:,:,nguard+1,lb) - &
        2.*solnData(isoln,:,:,nguard,lb)/delz**2
    endif

    if (nbrs(2) <= -20) then
      solnData(isrc,nguard+nxb,:,:,lb) = &
        solnData(isrc,nguard+nxb,:,:,lb) - &
        2.*solnData(isoln,nguard+nxb+1,:,:,lb)/delx**2
    endif

    if ((ndim >= 2) .and. (nbrs(4) <= -20)) then
      solnData(isrc,:,nguard+nyb,:,lb) = &
        solnData(isrc,:,nguard+nyb,:,lb) - &
        2.*solnData(isoln,:,nguard+nyb+1,:,lb)/dely**2
    endif

    if ((ndim == 3) .and. (nbrs(6) <= -20)) then
      solnData(isrc,:,:,nguard+nzb,lb) = &
        solnData(isrc,:,:,nguard+nzb,lb) - &
        2.*solnData(isoln,:,:,nguard+nzb+1,lb)/delz**2
    endif

  enddo

endif

!===============================================================================

return
end subroutine hg_init_src
