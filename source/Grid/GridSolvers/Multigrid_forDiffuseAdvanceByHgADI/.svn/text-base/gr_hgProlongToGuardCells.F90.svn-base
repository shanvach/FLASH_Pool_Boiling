!!****if* source/Grid/GridSolvers/Multigrid_forDiffuseAdvanceByHgADI/gr_hgProlongToGuardCells
!!
!! NAME
!!  gr_hgProlongToGuardCells
!!
!! SYNOPSIS
!!  
!!  call gr_hgProlongToGuardCells(integer(in) :: level,
!!                                integer(in) :: ifrom,
!!                                integer(in) :: ito,
!!                                integer(in) :: ichild)
!! DESCRIPTION
!!
!!  Prolong the boundary potentials present on a coarse level onto the
!!  face-guardcells of the children blocks at the next level.
!!  This occurs both on the faces and interiors of the parent blocks.
!!  For the faces, one interpolates the surface boundary values from the
!!  first layer of guardcells on the parent block.  For the second, one
!!  interpolates from the interior to get the approximated value at the
!!  bisection of the parent block.  These are put on the first layer of
!!  the child's boundary.  This scheme allows for seamless handling of
!!  external boundary values.  This routine is used in conjunction with a
!!  single-block solver to take an approximate solution on the coarse level
!!  and put it onto the fine levels.
!!  
!! ARGUMENTS
!!   
!!   ifrom  - the grid variable the data on the parent comes from
!!   ito    - the grid variable into which the face values are put
!!   level  - the level of the parents
!!   ichild - allows for only a single child block index to be filled
!!
!! NOTES
!!  
!!  Due to polynomial approximation's tendency to oscillate, the guardcell
!!  fills done here have two behaviors on the exterior blocks.  For
!!  the source term, the guardcells are extrapolated outwards.  For
!!  the residual terms, the guardcells are assumed 0 at the faces.
!! 
!!  This also increments the current level by calling gr_hgSetMaxLevel !
!!
!! SEE ALSO
!!  
!!  gr_hgProlongBndries
!!  gr_hgRestrict
!!
!!***

!!REORDER(5): unk

!==============================================================================

subroutine gr_hgProlongToGuardCells(level, ifrom, ito, ichild)

#include "Flash.h"
#include "Multigrid.h"
#include "constants.h"

  use gr_hgData, ONLY: hg_myPE, nbbuf_prolong, nmax1, nmax2, n1off, n2off, n3off,  &
       & Pns, Px, Py, Pz, Pud, Pew, send_prolong_req, send_prolong_data, recv_prolong_data, &
       & gr_hgSaveNodetype, gr_hgSolnIndex, gr_hgMeshRefineMax
  use Grid_interface, ONLY : Grid_fillGuardCells
  use Driver_interface, ONLY : Driver_abortFlash
  use paramesh_interfaces, ONLY : amr_prolong
  use tree, ONLY : lnblocks,parent,child,neigh,nfaces,nchild,nodetype,lrefine
  use physicaldata, ONLY : unk
  use workspace, ONLY : work

  use Timers_interface, ONLY : Timers_start, Timers_stop
  use gr_hgInterface, ONLY: gr_hgBndry, gr_hgSetMaxLevel,gr_hgLevelMultiplyScalar

  implicit none
  
#include "Flash_mpi.h"
  
  integer, intent(in)          :: ifrom, ito, level, ichild
  
  integer                      :: b, c, h, i, j, k, ierr1, ierr2, ierr3
  integer                      :: i1, i2, j1, j2, ii, jj, kk
  integer                      :: ierr, nsent
  integer                      :: blockID
  logical                      :: any_sent
  integer                      :: status(MPI_STATUS_SIZE)
  integer                      :: send_status(MPI_STATUS_SIZE,nchild*nbbuf_prolong)


!==============================================================================


  call Timers_start("gr_hgProlongToGuardCells")

  ! Use (EXCHANGE_WORK, CONTINUE_SERIES) under the assumption that we're
  ! calling gr_hgProlongToGuardCells immediately after a call to hg_solveLevel
  ! or hg_residual -- performance savings by not filling work() again
  if (ifrom == gr_hgSolnIndex) then
    call gr_hgBndry(level, ifrom, NGUARD, 0, MG_EXCHANGE_WORK, &
                    MG_STANDALONE, .true.)  ! do extrapolation across the boundary
  else
    call gr_hgBndry(level, ifrom, NGUARD, 0, MG_EXCHANGE_WORK, &
                    MG_STANDALONE, .false.)
  endif

  if (level < gr_hgMeshRefineMax) then
     call gr_hgSetMaxLevel(level+1)    ! This is a new FLASH3 routine

     call amr_prolong (hg_myPE, 2, NGUARD)
  
!     call gr_hgLevelMultiplyScalar(level, ifrom, 1.0, MG_NODES_ALL_NODES)
     if (ifrom == gr_hgSolnIndex) then
        call gr_hgBndry(level+1, ifrom, NGUARD, 0, MG_UPDATE_UNK, &
                    MG_END_SERIES, .true.)  ! do extrapolation across the boundary
     else
        call gr_hgBndry(level+1, ifrom, NGUARD, 0, MG_UPDATE_UNK, &
                    MG_END_SERIES, .false.)
     end if
!     call Grid_fillGuardCells( WORK, ALLDIR, &
!          doEos=.FALSE.,makeMaskConsistent=.FALSE.)
  end if

  call Timers_stop("gr_hgProlongToGuardCells")
  
  !===================================================================

  return
end subroutine gr_hgProlongToGuardCells
