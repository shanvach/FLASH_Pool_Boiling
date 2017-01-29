!!****if* source/Grid/GridSolvers/HYPRE_KPD/gr_hypreSetupSolver
!!
!!  NAME 
!!
!! gr_hypreSetupSolver
!!
!!  SYNOPSIS
!!
!!  call gr_hypreSetupSolver()
!!
!!
!!  DESCRIPTION 
!! This routine sets up the HYPRE solver object 
!! with it's associated PC. This routine is called only once.
!! 
!!
!! ARGUMENTS
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES:
!!  The solvers are present in HYPRE.We provide support for multiple solvers, however
!!  each of these solvers have too many knobs to tweak, with the usability in mind a
!!  few parameters have been selected which works fine for standard diffusion problems.
!!  
!!  However no way are these parameters the best for all problems, if needed additional
!!  HYPRE interfaces can be exposed in this routine.
!!
!!  For now PCG, BiCGSTAB, GMRES, AMG are exposed as solvers.
!!          AMG, ILU can be used as PC.
!!
!!  Again each of these solvers/pc have variety of parameters which are problem specific
!!  and need to be tweaked based on the problem being solved.
!!
!!-----------------------------------------------------------------------
!!     These settings will give consistent results across
!!     different distributions / domain decompositions.
!!     NOTE: Information obtained from HYPRE team and should be used
!!           only for debug purposes.
!!-----------------------------------------------------------------------     
!!     call HYPRE_BoomerAMGCreate(gr_hypreSolver, ierr)
!!     call HYPRE_BoomerAMGSetCoarsenType(gr_hypreSolver,7, ierr) ! CLJP coarsening
!!     call HYPRE_BoomerAMGSetRelaxType(gr_hypreSolver, 0, ierr) ! Jacobi-0, Chebyshev-16, FCF Jacobi-17, or l1 Jacobi-18
!!     call HYPRE_BoomerAMGSetMaxIter(gr_hypreSolver, gr_hypreMaxIter, ierr)
!!     call HYPRE_BoomerAMGSetTol(gr_hypreSolver, gr_hypreRelTol, ierr)  ! conv. tolerance     
!!
!!***

!!REORDER(4): solnVec

subroutine gr_hypreSetupSolver ()
  
  use gr_hypreData,  ONLY : gr_hypreSolver, gr_hyprePC,         &
                            gr_hypreRelTol, gr_hypreMaxIter,    &
                            gr_hyprePCType, gr_hypreSolverType, &    
                            gr_hyprePrintSolveInfo, gr_hypreInfoLevel, &
                            gr_hypreUse2Norm
  use Grid_data,     ONLY : gr_meshComm
  use Timers_interface, ONLY : Timers_start, Timers_stop  

  
  implicit none
  
#include "Flash.h"
#include "constants.h"   
#include "HYPREf.h"
  
  integer ::  precond_id, ierr

  integer :: sai_max_levels
  real    :: sai_threshold
  real    :: sai_filter     
  integer :: sai_sym

  call Timers_start("gr_hypreSetupSolver")     
  
  
  !!-----------------------------------------------------------------------
  !!     1.  Setup HYPRE Preconditioner object, if one is selected.
  !!-----------------------------------------------------------------------       

!!$  /*------------------------------------------------------------
!!$  * The precond_id flags mean :
!!$  *  0 - no preconditioner
!!$  *  1 - set up a ds preconditioner
!!$  *  2 - set up an amg preconditioner
!!$  *  3 - set up a pilut preconditioner
!!$  *  4 - set up a parasails preconditioner
!!$  *  5 - set up a euclid preconditioner (fix F90 interface in HYPRE).
!!$  *------------------------------------------------------------*/
  
  precond_id = 0
  
  select case (gr_hyprePCType)         
     
  case (HYPRE_AMG)     
     
     precond_id = 2     
     
      call HYPRE_BoomerAMGCreate(gr_hyprePC, ierr)       
      call HYPRE_BoomerAMGSetMaxIter(gr_hyprePC, 1, ierr)
      call HYPRE_BoomerAMGSetTol(gr_hyprePC, 0.0, ierr)    
     !call HYPRE_BoomerAMGSetCoarsenType(gr_hyprePC, 6, ierr)   ! Original... Bad 
!      call HYPRE_BoomerAMGSetCoarsenType(gr_hyprePC, 8, ierr)   ! KPD... Good
      call HYPRE_BoomerAMGSetRelaxType(gr_hyprePC, 6, ierr)     
      call HYPRE_BoomerAMGSetNumSweeps(gr_hyprePC, 1, ierr)      

     call HYPRE_BoomerAMGSetTruncFactor(gr_hyprePC, 0.5, ierr)

      !- kpd - Suggested From LLNL --
     call HYPRE_BoomerAMGSetCoarsenType(gr_hyprePC, 8, ierr)    ! 10 or 8... 6=BAD
     call HYPRE_BoomerAMGSetInterpType(gr_hyprePC, 6, ierr)     !  6 or 7  def=0
     call HYPRE_BoomerAMGSetPMaxElmts(gr_hyprePC, 4, ierr)      !  4 or 5  def=0
     call HYPRE_BoomerAMGSetAggNumLevels(gr_hyprePC, 1, ierr)   !  1 or 2  def=0


     !!! ADDITIONAL STUFF FROM ROB.
     !call HYPRE_BoomerAMGCreate(gr_hyprePC, ierr)         
     !call HYPRE_BoomerAMGSetMaxIter(gr_hyprePC, 1, ierr)
     !call HYPRE_BoomerAMGSetTol(gr_hyprePC, 0.0, ierr)    
     !call HYPRE_BoomerAMGSetRelaxType(gr_hyprePC, 6, ierr)     
     !call HYPRE_BoomerAMGSetNumSweeps(gr_hyprePC, 1, ierr)     
     !call HYPRE_BoomerAMGSetCoarsenType(gr_hyprePC, 10, ierr)
     !call HYPRE_BoomerAMGSetInterpType(gr_hyprePC, 6, ierr)
     !call HYPRE_BoomerAMGSetTruncFactor(gr_hyprePC, 0.5, ierr)
     !call HYPRE_BoomerAMGSetPMaxElmts(gr_hyprePC, 4, ierr)
     
!!$  !! For 3D
!!$  call HYPRE_ParCSRHybridSetPMaxElmts(amg_solver, 4, ierr)
     
     if (gr_hyprePrintSolveInfo) call HYPRE_BoomerAMGSetPrintLevel(gr_hyprePC,gr_hypreInfoLevel, ierr)           
     
!  case (HYPRE_PILUT)
!
!     precond_id = 3
!
!     print*,"Using the PILUT PreCond"
!     call HYPRE_ParCSRPilutCreate(gr_meshcomm, gr_hyprePC, ierr)
!!     call HYPRE_ParCSRPilutSetMaxIter(  )
!!     call HYPRE_

  case (HYPRE_PARASAILS)
     
     precond_id = 4
     
     sai_max_levels = 1
     !sai_threshold  = 0.1
     sai_threshold  = 0.05
     sai_filter     = 0.05
     sai_sym        = 1     
     
     call HYPRE_ParaSailsCreate(gr_meshcomm, gr_hyprePC, ierr)     
     call HYPRE_ParaSailsSetParams(gr_hyprePC, sai_threshold, sai_max_levels,ierr)
     call HYPRE_ParaSailsSetFilter(gr_hyprePC, sai_filter,ierr);
     call HYPRE_ParaSailsSetSym(gr_hyprePC, sai_sym,ierr)
     
     if (gr_hyprePrintSolveInfo) call HYPRE_ParaSailsSetLogging(gr_hyprePC,gr_hypreInfoLevel,ierr)         
     
  case (HYPRE_ILU)
     
     precond_id = 5 
     
     call HYPRE_EuclidCreate(gr_meshcomm, gr_hyprePC, ierr)             
     
     !! LEVELS OF FILL IN 
     call HYPRE_EuclidSetLevel(gr_hyprePC, 0, ierr)     

     !! 1 -> BILU or 0 -> PILU
     !call HYPRE_EuclidSetBJ(gr_hyprePC, 1, ierr)     
     call HYPRE_EuclidSetBJ(gr_hyprePC, 1, ierr)     

!     call HYPRE_EuclidSetStats(gr_hyprePC, 1, ierr)
     
     !if (gr_hyprePrintSolveInfo) call HYPRE_EuclidSetStats(gr_hyprePC,gr_hypreInfoLevel,ierr)
     
  end select
  
  !!-----------------------------------------------------------------------
  !!     2.  Setup HYPRE solver object, associate a PC with it.
  !!-----------------------------------------------------------------------    


  select case (gr_hypreSolverType)
     
  case (HYPRE_PCG)     

!print*,"Setting Up the PCG"

     call HYPRE_parCSRPCGCreate(gr_meshComm, gr_hypreSolver, ierr)     
     call HYPRE_ParCSRPCGSetTol(gr_hypreSolver, gr_hypreRelTol, ierr)          
     call HYPRE_ParCSRPCGSetMaxIter(gr_hypreSolver, gr_hypreMaxIter, ierr)       
     ! call HYPRE_ParCSRPCGSetTwoNorm(gr_hypreSolver, gr_hypreUse2Norm, ierr)
     call HYPRE_ParCSRPCGSetPrecond(gr_hypreSolver, precond_id, gr_hyprePC, ierr)

     if (gr_hyprePrintSolveInfo) call HYPRE_ParCSRPCGSetLogging(gr_hypreSolver, gr_hypreInfoLevel,ierr)     
     
  case (HYPRE_BICGSTAB)     
     call HYPRE_parCSRBiCGSTABCreate(gr_meshComm, gr_hypreSolver, ierr)     
     call HYPRE_ParCSRBiCGSTABSetTol(gr_hypreSolver, gr_hypreRelTol, ierr)         
     call HYPRE_ParCSRBiCGSTABSetMaxIter(gr_hypreSolver, gr_hypreMaxIter, ierr)  
     call HYPRE_ParCSRBiCGSTABSetPrecond(gr_hypreSolver, precond_id, gr_hyprePC, ierr)   

     if (gr_hyprePrintSolveInfo) call HYPRE_ParCSRBiCGSTABSetLogging(gr_hypreSolver,gr_hypreInfoLevel,ierr)
     
  case (HYPRE_AMG)     
     call HYPRE_BoomerAMGCreate(gr_hypreSolver, ierr)
!     call HYPRE_BoomerAMGSetCoarsenType(gr_hypreSolver, 6, ierr)            ! Falgout coarsening
     call HYPRE_BoomerAMGSetRelaxType(gr_hypreSolver, 3, ierr)              ! G-S/Jacobi hybrid relaxation
     call HYPRE_BoomerAMGSetNumSweeps(gr_hypreSolver, 1, ierr)              ! Sweeeps on each level
     call HYPRE_BoomerAMGSetMaxLevels(gr_hypreSolver, 20, ierr)             ! maximum number of levels */
     call HYPRE_BoomerAMGSetTol(gr_hypreSolver, gr_hypreRelTol, ierr)       ! conv. tolerance
     call HYPRE_BoomerAMGSetMaxIter(gr_hypreSolver, gr_hypreMaxIter, ierr)

      !- kpd - Suggested From LLNL --
     call HYPRE_BoomerAMGSetCoarsenType(gr_hypreSolver, 8, ierr)    ! 10 or 8... 6=BAD
     call HYPRE_BoomerAMGSetInterpType(gr_hypreSolver, 6, ierr)     !  6 or 7  def=0
     call HYPRE_BoomerAMGSetPMaxElmts(gr_hypreSolver, 4, ierr)      !  4 or 5  def=0
     call HYPRE_BoomerAMGSetAggNumLevels(gr_hypreSolver, 1, ierr)   !  1 or 2  def=0
     
     if (gr_hyprePrintSolveInfo) then
        call HYPRE_BoomerAMGSetPrintLevel(gr_hypreSolver,gr_hypreInfoLevel, ierr)      
        call HYPRE_BoomerAMGSetLogging(gr_hypreSolver,gr_hypreInfoLevel, ierr)    
     endif
     
  case (HYPRE_GMRES)          
     call HYPRE_ParCSRGMRESCreate(gr_meshComm, gr_hypreSolver, ierr)     
     call HYPRE_ParCSRGMRESSetKDim(gr_hypreSolver, 30, ierr)
     call HYPRE_ParCSRGMRESSetMaxIter(gr_hypreSolver, gr_hypreMaxIter, ierr)
     call HYPRE_ParCSRGMRESSetTol(gr_hypreSolver, gr_hypreRelTol, ierr)     
     call HYPRE_ParCSRGMRESSetPrecond(gr_hypreSolver, precond_id, gr_hyprePC, ierr)        

     if (gr_hyprePrintSolveInfo) call HYPRE_ParCSRGMRESSetLogging(gr_hypreSolver,gr_hypreInfoLevel,ierr)         
     
  case (HYPRE_SPLIT)     
     call HYPRE_SStructSplitCreate(gr_meshComm, gr_hypreSolver, ierr)
     call HYPRE_SStructSplitSetStructSolv(gr_hypreSolver, 17, ierr)
     call HYPRE_SStructSplitSetMaxIter(gr_hypreSolver, gr_hypreMaxIter, ierr)
     call HYPRE_SStructSplitSetTol(gr_hypreSolver, gr_hypreRelTol, ierr)     


  case (HYPRE_HYBRID)     
    
     call HYPRE_ParCSRHybridCreate(gr_hypreSolver, ierr)
     
     !! Set the convergence tolerance for the Krylov solver. The default is 1.e-7.
     call HYPRE_ParCSRHybridSetTol(gr_hypreSolver, gr_hypreRelTol, ierr)
     
     !! Set the maximal number of iterations for the diagonally preconditioned solver.
     call HYPRE_ParCSRHybridSetDSCGMaxIte(gr_hypreSolver,gr_hypreMaxIter, ierr) 
     
     !! Set the maximal number of iterations for the AMG preconditioned solver
     call HYPRE_ParCSRHybridSetPCGMaxIter(gr_hypreSolver,gr_hypreMaxIter, ierr) 
     
     !! Set the desired solver type. There are the following options:
     !! 1 PCG (default)
     !! 2 GMRES
     !! 3 BiCGSTAB     
     call HYPRE_ParCSRHybridSetSolverType(gr_hypreSolver, 1, ierr)
     

     if (gr_hyprePrintSolveInfo) then
        
        !! Set logging parameter (default: 0, no logging)
        call HYPRE_ParCSRHybridSetLogging(gr_hypreSolver, gr_hypreInfoLevel, ierr)

        !! Set print level (default: 0, no printing)
        call HYPRE_ParCSRHybridSetPrintLevel(gr_hypreSolver, gr_hypreInfoLevel, ierr)
        
     end if     
     
     !! (Optional) Defines a truncation factor for the interpolation. The default is 0.
     call HYPRE_ParCSRHybridSetTruncFacto(gr_hypreSolver, 0.5, ierr)
     
     !! (Optional) Defines the maximal number of elements per row for the interpolation. The default is 0.
     !! Does not have a F90 Interface.
!!$     call HYPRE_ParCSRHybridSetPMaxElmts(gr_hypreSolver, 4, ierr)
     
     !! (Optional) Defines which parallel coarsening algorithm is used.
     call HYPRE_ParCSRHybridSetCoarsenTyp(gr_hypreSolver,10,ierr)
     
     !! (Optional) Sets the number of sweeps. On the Finest level, the up and the down cycle the number of sweeps
     !! are set to num sweeps and on the coarsest level to 1. The default is 1.
     call HYPRE_ParCSRHybridSetNumSweeps (gr_hypreSolver, 1, ierr)     
     
     !! (Optional) Defines the smoother to be used. It uses the given smoother on the fine grid, the up and the
     !! down cycle and sets the solver on the coarsest level to Gaussian elimination (9). The default is Gauss-Seidel(3).
     call HYPRE_ParCSRHybridSetRelaxType(gr_hypreSolver, 6, ierr)
     
  end select 
  
  
  call Timers_stop("gr_hypreSetupSolver")    

  
  return
  
end subroutine gr_hypreSetupSolver
