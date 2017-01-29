!*******************************************************************************

!  Routine:     mg_initSlv()

!  Description: Multigrid initialization routine.


  subroutine gr_mgInitSlv(bndTypes)

!===============================================================================

#include "Flash.h"

  use mg_common

  use tree, only : nodetype,newchild,lrefine,maxblocks_tr

  use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                                Grid_getListOfBlocks 

  use Grid_data, ONLY : gr_geometry

  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none
#include "Multigrid.h"
#include "constants.h"
include "mpif.h"

  integer, intent(in) :: bndTypes(6)

  integer            :: lb, i, lmin, lmax, ierr, lnblocks2
  integer, parameter :: MAXDIM2 = 3
  integer            :: nbr_blks(2*MAXDIM2)

  logical, save :: FirstCall = .true.
  integer       :: igeom
  logical       :: bnd_is_valid

  real, pointer, dimension(:,:,:,:) :: unkt
  integer, pointer, dimension(:)  :: nodetype2
  logical, pointer, dimension(:)  :: newchild2

  integer, parameter :: nxb = NXB
  integer, parameter :: nyb = NYB
  integer, parameter :: nzb = NZB
  integer, parameter :: nguard = NGUARD
  integer, parameter :: ndim = NDIM
  integer, parameter :: k2d = K2D
  integer, parameter :: k3d = K3D

  integer :: eachboundary


!===============================================================================


! INITILIAZE PARAMESH multigrid specific data
call amr_mg_init()

call Grid_getLocalNumBlks(lnblocks2)

do lb = 1, lnblocks2
   nodetype_save(lb) = nodetype(lb)
   newchild_save(lb) = newchild(lb)
end do

! Determine the smallest and largest levels of refinement in the mesh.
lmax = -1
lmin = 10000000
do i = 1, lnblocks2
  if (nodetype(i) == 1) then
    lmin = min( lrefine(i), lmin )
    lmax = max( lrefine(i), lmax )
  endif
enddo

call mpi_allreduce ( lmin, mesh_lrefmin, 1, MPI_INTEGER, MPI_MIN, & 
                     MPI_COMM_WORLD, ierr )
call mpi_allreduce ( lmax, mesh_lrefmax, 1, MPI_INTEGER, MPI_MAX, & 
                     MPI_COMM_WORLD, ierr )


! Assign Boundary condition types to gr_mgBndTypes, same as hg solver:
  do eachBoundary = 1, 2*NDIM

     gr_mgBndTypes(eachBoundary) = bndTypes(eachBoundary)

!!$     if (bndTypes(eachBoundary) == gr_mgbcTypes(eachBoundary) .OR. suppressPfft .OR. &
!!$        (bndTypes(eachBoundary) == MG_BND_GIVENVAL .AND.                             &
!!$         gr_mgbcTypes(eachBoundary)==MG_BND_DIRICHLET) ) then
!!$        gr_mgBndTypes(eachBoundary) = bndTypes(eachBoundary)
!!$     else
!!$        if (gr_myPe .eq. 0) then
!!$           write(*,*) 'gr_mgSolve Error: Boundary Conditions for Poisson Solver is inconsistent.'
!!$           write(*,*) 'gr_mgSolve Error: direction=',eachBoundary
!!$           write(*,*) 'gr_mgSolve Error: gr_mgbcTypes(direction) =',gr_mgbcTypes(eachBoundary)
!!$           write(*,*) 'gr_mgSolve Error: bndTypes(direction)     =',bndTypes(eachBoundary)
!!$        endif
!!$        call Driver_abortFlash('gr_mgSolve Error: BC type argument inconsistent')
!!$     end if
  end do

! Assign value to mg_bnd_cond: case 0, substract mean from source, 
!                              case 1, subtract given value of solution as
!                                      in Dirichlet BCs.

  mg_bnd_cond = 0 !Assume all BCs are Periodic or Neumann
  do  eachBoundary = 1, 2*NDIM
    if ((bndTypes(eachBoundary)==MG_BND_GIVENVAL) .or. &
        (bndTypes(eachBoundary)==MG_BND_DIRICHLET)   ) &
        mg_bnd_cond = 1 ! Case Some BC is DIRICHLET-GIVENVAL
  enddo

!===============================================================================

return
end
