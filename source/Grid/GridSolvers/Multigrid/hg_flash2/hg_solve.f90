!******************************************************************************

! Routine:      hg_solve

! Description:  Solve the Poisson equation on a FLASH/PARAMESH mesh using a
!               variant of the method of Huang & Greengard (2000, SIAM J. Sci.
!               Comput., 21, 1551).  The originally published method has been
!               modified to work with the PARAMESH mesh structure and zone-
!               averaged rather than corner values.

! Parameters:   isrc        Variable index for source function (must be valid
!                           on all blocks)
!               isoln       Variable index for solution
!               isls        Variable index for single-layer strength (residual)
!               icorr       Variable index for correction
!               SolveBlock  A routine to call a single-block direct Poisson
!                           solver once boundary values have been set.
!               bnd_type    Boundary condition type (types defined in
!                           mg_common)
!               src_fact    Factor multiplying source term on RHS of equation
!                           being solved


subroutine hg_solve(isrc, isoln, isls, icorr, SolveBlock, bnd_type, src_fact)

!==============================================================================

use dBase, ONLY: dBasePropertyInteger, nguard
use mg_common
use runtime_parameters
use perfmon

implicit none

integer, intent(in) :: isrc, isoln, isls, icorr, bnd_type, src_fact
external               SolveBlock

integer             :: n, m, ierr
real                :: norm_lhs, norm_rhs, norm, norm_old

integer, save       :: max_mg_corrections, MyPE, MasterPE
real, save          :: max_mg_residual_norm
logical, save       :: first_call = .true., mgrid_print_norm

!==============================================================================

call timer_start("hg_solve")

! Generic initializations.

if (first_call) then
  MyPE     = dBasePropertyInteger("MyProcessor")
  MasterPE = dBasePropertyInteger("MasterProcessor")
  call get_parm_from_context("max_mg_corrections", max_mg_corrections)
  call get_parm_from_context("max_mg_residual_norm", max_mg_residual_norm)
  call get_parm_from_context("mgrid_print_norm", mgrid_print_norm)
  first_call = .false.
endif

! Initialize source data.

mg_soln_index = isoln
mg_bnd_cond = bnd_type
call hg_init_src(isrc, isoln)
call mg_norm(0, isrc, norm_rhs, MG_NODES_LEAF_ONLY)

! Coarse grid interpolation step.  Solve on coarse grid, interpolate boundary
! conditions for next finer level, solve on that level, etc.

do m = 1, mesh_lrefmax
  call hg_set_ext_bndries(m, isoln, 0)
  call hg_solve_level(m, isrc, isoln, SolveBlock, MG_NODES_ALL_NODES)
  call hg_residual(m, isrc, isoln, isls)
  call hg_prolong_bndries(m, isoln, isoln, 0)
enddo

! Correction step.  Restrict residuals from finer levels to coarser levels.
! Solve for correction on these levels and interpolate boundary conditions to
! finer levels.  Solve for corrections there and apply.  Repeat.  Repeat these
! correction steps until the desired residual norm is achieved.

do n = 0, max_mg_corrections

  call mg_norm(0, isls, norm_lhs, MG_NODES_LEAF_ONLY)
  if ((mgrid_print_norm) .and. (MyPE == MasterPE)) then
    if (n == 0) then
      write(*,'(A,I3,A,ES12.5)') &
        'hg_solve: iter ', n, ': norm(res)/norm(src) = ', &
        norm_lhs/norm_rhs
    else
      write(*,'(A,I3,2(A,ES12.5))') &
        'hg_solve: iter ', n, ': norm(res)/norm(src) = ', &
        norm_lhs/norm_rhs, ' convg fact = ', norm_lhs/norm_old
    endif
  endif
  norm_old = norm_lhs

  if (norm_lhs/norm_rhs <= max_mg_residual_norm) exit

  do m = mesh_lrefmax-1, 1, -1
    call hg_restrict(m+1, isls, isls)
  enddo

  call hg_set_ext_bndries(1, icorr, 0)
  call hg_solve_level(1, isls, icorr, SolveBlock, MG_NODES_ALL_NODES)
  call hg_level_add(1, isoln, icorr, MG_NODES_LEAF_ONLY)
  call hg_prolong_bndries(1, icorr, icorr, 0)

  do m = 2, mesh_lrefmax
    call hg_set_ext_bndries(m, icorr, 0)
    call hg_solve_level(m, isls, icorr, SolveBlock, MG_NODES_ALL_NODES)
    call hg_residual(m, isls, icorr, isls)
    call hg_level_add(m, isoln, icorr, MG_NODES_LEAF_ONLY)
    call hg_prolong_bndries(m, icorr, icorr, 0)
  enddo

enddo

! Multiply solution by the source term factor.

do m = 1, mesh_lrefmax
  call hg_level_smultiply(m, isoln, src_fact, MG_NODES_LEAF_ONLY)
enddo

! Add back source average to source function (to compensate for its subtraction
! when using periodic/Neumann boundary conditions).

do m = 1, mesh_lrefmax
  call hg_level_sadd(m, isrc, src_avg, MG_NODES_LEAF_ONLY)
enddo

! Leave boundary zones properly updated.

call mg_bndry (0, isoln, nguard, MG_NODES_LEAF_ONLY, UPDATE_UNK, &
               STANDALONE, 1)

! Make sure we leave PARAMESH in a sane state.

!call mg_restore_nodetypes()

call timer_stop("hg_solve")

!==============================================================================

return
end subroutine hg_solve
