!*******************************************************************************

!  Routine:     multigrid()

!  Description: Main multigrid Poisson solver.  The solution is preconditioned
!               on the first call; afterward, the existing solution is taken
!               as the initial guess.

!               multigrid multiplies unk(isrc_dens) by poisfact and puts the
!               result into unk(img_src).  If the img_bc arg is 0 (dirichlet/
!               given-value), take boundary values for potential from boundaries
!               of unk(1st arg).  Then solve.  On output, solution goes into
!               unk(1st arg); unk(2nd arg) is not modified (unless 2nd
!               arg=img_src).

!  Parameters:  isrc_dens       Index for source array.  This is taken to be
!                               the density field; the source array to be used
!                               as the right-hand side of the Poisson equation
!                               is computed from this.
!               img_soln        Index for solution array.  The solution is
!                               written directly into this variable.
!               poisfact        Constant Poisson factor.
!               img_src, img_res, img_corr, img_temp
!                               Indices for work arrays to hold, respectively,
!                               the source (right-hand side), residual,
!                               correction, and temporary values.
!               img_bc          Boundary condition to apply to all boundaries.
!                                 0 = Isolated boundaries using James' algorithm
!                                 1 = Periodic boundaries
!                                 2 = Dirichlet boundaries
!                                 3 = Neumann boundaries


subroutine multigridHYPRE (isrc_dens, img_soln, poisfact, img_src, img_res, & 
                      img_corr, img_temp, img_temp2, bc_types,mg_solve,  &
                      mg_residual, mg_residualMG, mg_relax_RBGS, &
                      mg_relax_HYPRE, &
                      img_denx, mg_relax)

!===============================================================================

#include "Flash.h"

  use mg_common, ONLY: ile, iue, jle, jue, kle, kue, &
     mg_bnd_cond, mesh_lrefmax, interp_work,interp_mask_work_mg, &
     interp_mask_work_save,nodetype_save, &
     newchild_save, gr_mgDiffOpDiscretize, mesh_lrefmin

  use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                                Grid_getListOfBlocks, &
                                Grid_getBlkPtr,       &
                                Grid_releaseBlkPtr

  use Grid_data, ONLY : gr_meshMe

  use gr_mgInterface, ONLY : gr_mgInitSlv,gr_mgInitSrc, &
    gr_mgNorm, mg_restore_nodetypes, gr_mgBndry

  use gr_interface, ONLY : gr_findMean

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  
  use workspace, ONLY: work, interp_mask_work
  use paramesh_dimensions, only : nguard_work,nvar_work

  use tree, only : maxblocks_tr

implicit none
#include "Multigrid.h"
#include "constants.h"

integer       :: isrc_dens, img_soln, img_src, img_res, img_corr, img_temp, &
                 img_temp2, img_denx
integer       :: bc_types(6)
real          :: poisfact

external mg_solve, mg_residual, mg_residualMG, mg_relax_RBGS, &
         mg_relax_HYPRE, mg_relax

integer :: blockCount
integer :: blockList(MAXBLOCKS)

logical       :: done
integer       :: i, j, k, ii, jj, kk, lb, lnblocks, level
real          :: res_norm_ratio, res_norm_change, norm_old, norm_new, norm_rhs
real, save    :: mgrid_max_residual_norm, mgrid_max_iter_change
integer, save :: mgrid_max_vcycles, MyPE, MasterPE
logical, save :: mgrid_print_norm
logical, save :: first_call = .true.
logical, save :: first_step = .true.

integer :: blockID,img_bc

real, pointer, dimension(:,:,:,:) :: unk

real :: mean_soln

!- kpd - Added to check on variables that are passed into the routines...
real, pointer, dimension(:,:,:,:) :: solnData

!===============================================================================

! Initialize
!call timer_start("multigrid")

call Grid_getListOfBlocks(LEAF,blockList,blockCount)

!-----------------------------------------------------------------------------------
!- kpd - This is only called during the first time step CURRENTLY
!-----------------------------------------------------------------------------------
if (first_call) then

  call RuntimeParameters_get('mgrid_max_residual_norm',    mgrid_max_residual_norm)
  call RuntimeParameters_get('mgrid_max_iter_change',    mgrid_max_iter_change)
  call RuntimeParameters_get('mgrid_max_vcycles',    mgrid_max_vcycles)
  call RuntimeParameters_get('mgrid_print_norm',    mgrid_print_norm)

  MyPE = gr_meshMe 
  MasterPE = MASTER_PE

  allocate(interp_mask_work_mg(nvar_work))

  ! This statement is missing in the non ppft version of multigrid !!!!
  call RuntimeParameters_get("gr_mgDiffOpDiscretize",gr_mgDiffOpDiscretize)
  !    -------

  !- kpd - The discretization scheme is set to 2nd or 4th order central
  !--------------------------------------------------------------------
  select case(gr_mgDiffOpDiscretize)
  case(2) ! 2nd Order Central differences
    interp_work = 2
  case(4) ! 4th Order Central differences
    interp_work = 4
    write(*,*) 'IN Multigrid START gr_mgDiffOpDiscretize = 4'
  end select

  interp_mask_work_mg = interp_work
  allocate(interp_mask_work_save(nvar_work))

  allocate(nodetype_save(maxblocks_tr))
  allocate(newchild_save(maxblocks_tr))

  first_call = .false.
endif
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------

 interp_mask_work_save = interp_mask_work;
 interp_mask_work = interp_mask_work_mg;

!-------------------------------------
!- kpd - Initialize ParaMesh Multigrid
!-------------------------------------
call gr_mgInitSlv(bc_types)

!---------------------------------------------------------------------------------
!- kpd - Initialize the RHS of the Poisson Eqn (poisfact is 1.0), so 
!           img_src = poisfact * isrc_dens = DUST_VAR. 
!        img_soln is not changed (it is used fror BC setup)
!---------------------------------------------------------------------------------
call gr_mgInitSrc (isrc_dens, poisfact, img_src, img_soln)

!----------------------------------------------------------------------------------
!- kpd - Pressure (or DELP) is set to 0. This is only done during the 1st iteration
!----------------------------------------------------------------------------------
if (first_step) then
   call Grid_getListOfBlocks(LEAF,blockList,blockCount)                            !- kpd - This is only done for level 0 (or4).
   do lb = 1, blockCount                                                           !           NOT the level 5 blocks!
      blockID = blockList(lb)
      !print*,"KPD Initial Leaf Blocks Before mgCycle",blockID
      ! Point to blocks center vars:
      call Grid_getBlkPtr(blockID,unk,CENTER)

      !-------------------------------------------------------------------
      !- kpd - Initialize the solution (Pressure Error Correction) to Zero
      !-------------------------------------------------------------------
      unk(img_soln,:,:,:) = 0.

      ! Release pointers:
      call Grid_releaseBlkPtr(blockID,unk,CENTER)
   enddo
   first_step = .false.
endif


!---------------------------------------------------------------------------------
! - kpd - Iterate for the pressure solution. Convergence is checked at end of loop
!---------------------------------------------------------------------------------
done = .false.
i    = 0

!-------------------------------------------------------------------------------
!- kpd - Calculate the L2 Norm of the RHS (source), used for convergence testing
!        img_src is unchanged (intent - in)
!-------------------------------------------------------------------------------
call gr_mgNorm (0, img_src, norm_rhs, 1)                                          !- kpd - This is only done for level 0 (or 3).
               !^!                    ^                                           !           NOT the level 4 blocks!
               !Level #0              Leaf_Only
!---------------------------------------------------------------------
!- kpd - Calculate the Residual on this mg level. The residual will
!          will be used as the RHS for Poisson solution based on the
!          Residual Correction Method (Ae=r where r=b-Ap).
!        Note: The initial residual at the top of the V-cycle
!                 is calculated using the pressure from the 
!                 previous time step (img_soln), so the initial 
!                 residual is NOT the RHS (A*0=b, r=b)
!              The residual (img_res) is always zerod out heading
!                 into the residual call.
!---------------------------------------------------------------------
!- kpd - If a level is lower than lrefmin, than there are no blocks 
!          called in mg_residual (no harm, but wasted time, it should
!          go from mesh_lrefmax, mesh_lrefmin, -1)
!---------------------------------------------------------------------
!print*,"Before First of two mg_residual call in multigrid.F90 at level =",mesh_lrefmin,mesh_lrefmax
!!!do level = mesh_lrefmax, 1, -1


!- kpd - Assign mixture density to all leaf blocks
!!!!if (gr_meshMe .eq. 0) print*,"Using gr_mgAssignFden in multigridHYPRE"
!!!!call gr_mgAssignFden (mesh_lrefmax, img_denx, 0)

do level = mesh_lrefmax, mesh_lrefmin, -1

   !- kpd - img_res receives the updated soln for leaf blocks ONLY at all levels
   !                                              ==== ====== ====
   !                            RHS      LHS       RES
   !!!call mg_residual   (level, img_src, img_soln, img_res, 1, 1)
   !!!if (gr_meshMe .eq. 0) print*,"Using mg_residualMG instead of mg_residual in multigridHYPRE"
   call mg_residualMG (level, img_src, img_soln, img_res, 1, 1, img_denx, mesh_lrefmax)   !- kpd - ????


   !- kpd - It is IMPORTANT to do this b/c the nodetypes get altered in mg_residual
   call mg_restore_nodetypes (level)
enddo
!print*,"After First of two mg_residual call in multigrid.F90 at level =",mesh_lrefmin,mesh_lrefmax

!-----------------------------------------------------------------
!- kpd - Calculate the L2 Norm of the residual at the 0-leaf level
!           that was just calculated above in mg_residual
!        img_res is unchanged (intent - in)
!        Store norm as old norm headed into Poisson solve
!-----------------------------------------------------------------
call gr_mgNorm (0, img_res, norm_old, 1)
               !^!                    ^                                            !        
               !Level #0              Leaf_Only

!----------------------------------------------
!- kpd - Screen Output To Display AMR Levels...
!----------------------------------------------
if (gr_meshMe .eq. 0) then
   print*,"KPD - Entering V-Cycles with Min/Max no. of levels: ",mesh_lrefmin,mesh_lrefmax
end if

!-----------------------------------------------
! if the solution is not good enough, then cycle
!-----------------------------------------------
if (norm_old > mgrid_max_residual_norm * norm_rhs) then
   
   !call timer_start("mg_cycle")

   !----------------------------------------
   !- kpd - Iterate until convergence is met
   !----------------------------------------
   do while (.not. done)
  
      !---------------------------------------------------------------------------- 
      !- kpd - Perform one multigrid cycle beginning at a specified mesh level.
      !        *** The actual Poisson solver (mg_relax) is called in gr_mgCycle ***
      !        The updated pressure (img_soln) is the incoming pressure from the
      !           previous time step plus an error correction (img_corr)
      !        The error correction is determined from solving Ax=b with the 
      !           residual as the RHS(b)
      !        The initial error correction is the value from the previous time
      !           step, but this doesn't matter, b/c it is overwritten before
      !           it is used
      !        The RHS (img_src) and the residual (img_res) are unchanged.
      !---------------------------------------------------------------------------- 
      !print*,"KPD Entering gr_mgCycle in multigrid.F90",gr_meshMe
      call gr_mgCycleHYPRE (mesh_lrefmax,mesh_lrefmin, img_soln, img_src, img_res, img_corr, &
           img_temp, img_temp2, mg_solve, mg_residual, mg_residualMG, mg_relax_RBGS, &
           mg_relax_HYPRE, img_denx,mg_relax)
!      call gr_mgCycle (mesh_lrefmax, img_soln, img_src, img_res, img_corr, &
!           img_temp, img_temp2, mg_solve, mg_residual, mg_relax)
      !print*,"KPD Leaving gr_mgCycle in multigrid.F90",gr_meshMe

!********************************************************************************************
!********************************************************************************************
!********************************************************************************************

      !-------------------------------------------------------------------
      !- kpd - Compute the residual from the multigrid solve in gr_mgCycle
      !        This residual call is used for convergence check, or
      !           as the RHS if another Multigrid cycle is needed.
      !-------------------------------------------------------------------
      !!!do level = mesh_lrefmax, 1, -1
      do level = mesh_lrefmax, mesh_lrefmin, -1

         !!!call mg_residual   (level, img_src, img_soln, img_res, 1, 1)
         !!!if (gr_meshMe .eq. 0) print*,"Using mg_residualMG instead of mg_residual in multigridHYPRE"
         call mg_residualMG (level, img_src, img_soln, img_res, 1, 1, img_denx, mesh_lrefmax)   !- kpd - ????

         !- kpd - It is IMPORTANT to do this b/c the nodetypes get altered in mg_residual
         call mg_restore_nodetypes (level)

      enddo


      !------------------------------------------------------
      !- kpd - Compute the updated (new) norm of the residual
      !------------------------------------------------------
      call gr_mgNorm (0, img_res, norm_new, 1)

 
      res_norm_ratio  = norm_new / norm_rhs
      res_norm_change = (norm_new - norm_old) / norm_old
      
      if ((i > 0) .and. (res_norm_change > 0. .and. &
      norm_new > mgrid_max_residual_norm) .and. (MyPE == MasterPE)) &
           print *, "multigrid:  WARNING:  V-cycles not converging"

      i = i + 1

      !-----------------------------------------------------
      !- kpd - Output multigrid solver results to the screen
      !-----------------------------------------------------
      if ((mgrid_print_norm) .and. (MyPE == MasterPE)) then
      print*,"-----------------------------------------------------------------------------------------------------"
               if (norm_rhs > 0.e0) then
                  write(*,'(a,i4,3(a,es9.2))')                 &
                    ' Multigrid cycle # ', i,                  &
                    ' : res_norm_ratio = ', res_norm_ratio,    &
                    ' res norm = ',         norm_new,          &
                    ' new to old ratio = ', norm_new/norm_old
               else
                  write(*,'(a,i4,a,es9.2)')                    &
                     ' Multigrid cycle # ', i,                 &
                     ' : res_norm_change = ', res_norm_change
               end if
      print*,"-----------------------------------------------------------------------------------------------------"
      endif
      !print *, 'cycle ', i, ':  res_norm_ratio = ', res_norm_ratio

      !-----------------------------------------------------------
      !- kpd - Check convergence criterion, and exit if met (true)
      !-----------------------------------------------------------
      if ((i == mgrid_max_vcycles) .or. & 
           (norm_new <= mgrid_max_residual_norm * norm_rhs) .or. &
           (abs(res_norm_change) < mgrid_max_iter_change) .and.  &
           (i .gt. 0)) &
           done = .true.

      norm_old = norm_new

   enddo
   !---------------------------------------------------------------
   !---------------------------------------------------------------

   !call timer_stop("mg_cycle")

!********************************************************************************************
!********************************************************************************************
!********************************************************************************************

   !print*,"KPD - Entering mg_restore_nodetypes Before gr_mgBndry"

   !-------------------------------------------------------------
   ! Make sure the boundary zones are set properly
   !- kpd - Also, copy the pressure solution into the work array.
   !        I think this also has to do with AMR matching, etc.
   !-------------------------------------------------------------
   !do level = 1, mesh_lrefmax
   do level = mesh_lrefmin, mesh_lrefmax

      !- kpd - I put nodetype restoration in b/c paramesh was crashing for AMR with >1 level
      call mg_restore_nodetypes (level)

      call gr_mgBndry (i, img_soln, nguard_work, 1, MG_COPY_UNK_TO_WORK, MG_CONTINUE_SERIES)
   enddo

   ! call Grid_getLocalNumBlks(lnblocks)
   call Grid_getListOfBlocks(LEAF,blockList,blockCount)

   !print*,"KPD - Leaving mg_restore_nodetypes Before gr_mgBndry"

   !=====================================================================
   !- kpd - Here the final Pressure solution is stored in unk(img_soln, )
   !=====================================================================
   do lb = 1, blockCount
      blockID = blockList(lb)
      !print*,"KPD Final Leaf Blocks After mgCycle",blockID
      ! Point to blocks center vars:
      call Grid_getBlkPtr(blockID,unk,CENTER)
      do k = kle, kue
         do j = jle, jue
            do i = ile, iue

               unk(img_soln,i,j,k) = work(i,j,k,blockID,1)

               !print*,"HYPREsoln",lb,blockID,i,j,unk(img_soln,i,j,k),img_soln

            enddo
         enddo
      enddo
      ! Release pointers:
      call Grid_releaseBlkPtr(blockID,unk,CENTER)
   enddo


end if

 interp_mask_work = interp_mask_work_save

!print*,"KPD - Leaving multigridHYPRE"

!call timer_stop("multigrid")
!===============================================================================

return
end

