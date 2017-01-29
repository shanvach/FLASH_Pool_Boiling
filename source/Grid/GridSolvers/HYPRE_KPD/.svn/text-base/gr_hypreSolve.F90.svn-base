!!****if* source/Grid/GridSolvers/HYPRE_KPD/gr_hypreSolve
!!
!!  NAME 
!!
!!  gr_hypreSolve
!!
!!  SYNOPSIS
!!
!!  call gr_hypreSolve ()
!!
!!  DESCRIPTION 
!!      This routine solves AX=B using HYPRE
!!
!! ARGUMENTS
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!!
!!      Requires HYPRE library.
!!
!!***

!!REORDER(4): solnVec


subroutine gr_hypreSolve()
  
  use gr_hypreData,     ONLY : gr_hypreMatA, gr_hypreVecB, gr_hypreVecX, &
                               gr_hypreSolver, gr_hypreSolverType, &
                               gr_hyprePrintSolveInfo,gr_hypreMaxIter
  use Timers_interface, ONLY : Timers_start, Timers_stop  
  use Grid_data,        ONLY : gr_meshMe

  use Driver_data,      ONLY : dr_nstep

  use tree,             ONLY : grid_changed
  
  implicit none
  
#include "Flash.h"  
#include "constants.h"
#include "HYPREf.h"   

  integer, parameter :: range=SELECTED_INT_KIND(16)

  integer (KIND=range)::  parA
  integer (KIND=range)::  parb
  integer (KIND=range)::  parx 
  integer :: ierr

  integer :: num_iterations
  real    :: final_res_norm 

  character(len=32) :: matfile
  
  call Timers_start("gr_hypreSolve")     
  
  call HYPRE_SStructVectorAssemble(gr_hypreVecB, ierr)  
  
  num_iterations   = 0  
  final_res_norm   = 0.0
  
  if (gr_hypreSolverType == HYPRE_SPLIT) then     
     call HYPRE_SStructSplitSetup(gr_hypreSolver, gr_hypreMatA, gr_hypreVecB, gr_hypreVecX, ierr)
     call HYPRE_SStructSplitSolve(gr_hypreSolver, gr_hypreMatA, gr_hypreVecB, gr_hypreVecX, ierr)
  else  

     if (dr_nstep .eq. 1) then

        !matfile = 'amatb4slv'
        !matfile(10:10) = char(0)
        !call HYPRE_SStructMatrixPrint(matfile, gr_hypreMatA, 0, ierr)
     
        !matfile = 'bvecb4slv'
        !matfile(10:10) = char(0)
        !call HYPRE_SStructVectorPrint(matfile, gr_hypreVecB, 0, ierr)

     end if

     !!-----------------------------------------------------------------------
     !! Because we are using a PARCSR solver, we need to get the object
     !! of the matrix and vectors to pass in to the ParCSR solvers.
     !!-----------------------------------------------------------------------  
     call HYPRE_SStructMatrixGetObject(gr_hypreMatA, parA, ierr)
     call HYPRE_SStructVectorGetObject(gr_hypreVecB, parb, ierr)
     call HYPRE_SStructVectorGetObject(gr_hypreVecX, parx, ierr)
     
if (gr_meshMe .eq. 0) print*,"KPD - Solver Type",gr_hypreSolverType,HYPRE_PCG

     if (gr_hypreSolverType == HYPRE_PCG) then            
        
!        if (grid_changed .eq. 1 .OR. MOD(dr_nstep,200) .eq.0 ) then
!           if (gr_meshMe == 0) print*,"KPD - Grid Has Changed... HYPRE PCG Setup"
        if (gr_meshMe == 0) print*,"KPD - Entering HYPRE PCG Setup"
        call Timers_start("HYPRE_ParCSRPCGSetup")    
        call HYPRE_ParCSRPCGSetup(gr_hypreSolver, parA, parb, parx, ierr) 
        call Timers_stop("HYPRE_ParCSRPCGSetup")   
!        else
!           if (gr_meshMe == 0) print*,"KPD - Ignoring HYPRE PCG Setup, MOD=",MOD(dr_nstep,200)
!        end if        

        if (gr_meshMe == 0) print*,"KPD - Entering HYPRE PCG Solve"
        call Timers_start("HYPRE_ParCSRPCGSolve")    
        call HYPRE_ParCSRPCGSolve(gr_hypreSolver, parA, parb, parx, ierr)     
        call Timers_stop("HYPRE_ParCSRPCGSolve")        

        
     else  if (gr_hypreSolverType == HYPRE_BICGSTAB) then
        
        call Timers_start("HYPRE_ParCSRBiCGSTABSetup")    
        call HYPRE_ParCSRBiCGSTABSetup(gr_hypreSolver, parA, parb, parx, ierr)      
        call Timers_stop ("HYPRE_ParCSRBiCGSTABSetup")    
        
        call Timers_start("HYPRE_ParCSRBiCGSTABSolve")    
        call HYPRE_ParCSRBiCGSTABSolve(gr_hypreSolver, parA, parb, parx, ierr)     
        call Timers_stop("HYPRE_ParCSRBiCGSTABSolve")    
        
     else if(gr_hypreSolverType == HYPRE_AMG) then         
        
        call Timers_start("HYPRE_BoomerAMGSetup")
        call HYPRE_BoomerAMGSetup(gr_hypreSolver, parA, parb, parx, ierr)
        call Timers_stop("HYPRE_BoomerAMGSetup")
        
        call Timers_start("HYPRE_BoomerAMGSolve")
        call HYPRE_BoomerAMGSolve(gr_hypreSolver, parA, parb, parx, ierr)
        call Timers_stop("HYPRE_BoomerAMGSolve")                            
        
     else if(gr_hypreSolverType == HYPRE_GMRES) then   
        
        call Timers_start("HYPRE_ParCSRGMRESSetup")    
        call HYPRE_ParCSRGMRESSetup(gr_hypreSolver, parA, parb, parx, ierr)      
        call Timers_stop("HYPRE_ParCSRGMRESSetup")    
        
        call Timers_start("HYPRE_ParCSRGMRESSolve")    
        call HYPRE_ParCSRGMRESSolve(gr_hypreSolver, parA, parb, parx, ierr)
        call Timers_stop("HYPRE_ParCSRGMRESSolve")    
        
     else if(gr_hypreSolverType == HYPRE_HYBRID) then           
        
        call Timers_start("HYPRE_ParCSRHybridsetup")    
        call HYPRE_ParCSRHybridsetup(gr_hypreSolver,parA, parb, parx, ierr)
        call Timers_stop("HYPRE_ParCSRHybridsetup")    
        
        call Timers_start("HYPRE_ParCSRHybridSolve")    
        call HYPRE_ParCSRHybridSolve(gr_hypreSolver, parA, parb, parx, ierr)
        call Timers_stop("HYPRE_ParCSRHybridSolve")    
     end if
     
  end if

  if (gr_meshMe == 0) then
     
     select case (gr_hypreSolverType)
        
     case (HYPRE_PCG)        
200	FORMAT (ES11.4)
        call HYPRE_ParCSRPCGGetNumIterations(gr_hypreSolver, num_iterations, ierr)        
        call HYPRE_ParCSRPCGGetFinalRelative(gr_hypreSolver, final_res_norm, ierr)
        print*, "HYPRE PCG Num Iterations = ", num_iterations 
        print*, "HYPRE PCG Relative Residual Norm = ", final_res_norm
        WRITE(*,200) final_res_norm
        
     case (HYPRE_AMG)        
        call HYPRE_BoomerAMGGetNumIterations(gr_hypreSolver, num_iterations, ierr)
        call HYPRE_BoomerAMGGetFinalReltvRes(gr_hypreSolver, final_res_norm, ierr)
!!$        print*, "HYPRE AMG Num Iterations = ", num_iterations 
!!$        print*, "HYPRE AMG Relative Residual Norm = ", final_res_norm        
        
     case (HYPRE_BICGSTAB)                
        call HYPRE_ParCSRBICGSTABGetNumIter(gr_hypreSolver, num_iterations, ierr)
        call HYPRE_ParCSRBICGSTABGetFinalRel(gr_hypreSolver, final_res_norm, ierr)        
        print*, "HYPRE BiCGSTAB Num Iterations = ", num_iterations 
        print*, "HYPRE BiCGSTAB Relative Residual Norm = ", final_res_norm
        
     case (HYPRE_GMRES)                
        call HYPRE_ParCSRGMRESGetNumIteratio(gr_hypreSolver, num_iterations, ierr)
        call HYPRE_ParCSRGMRESGetFinalRelati(gr_hypreSolver, final_res_norm, ierr)        
        print*, "HYPRE GMRES Num Iterations = ", num_iterations 
        print*, "HYPRE GMRES Relative Residual Norm = ", final_res_norm
        
     case (HYPRE_HYBRID)           
        
        call HYPRE_ParCSRHybridGetNumIterati(gr_hypreSolver, num_iterations, ierr)
!!$        print*, "HYPRE HYBRID Num Iterations      = ", num_iterations 
        
        call HYPRE_ParCSRHybridGetPCGNumIter (gr_hypreSolver, num_iterations, ierr)
!!$        print*, "HYPRE HYBRID (AMG)Num Iterations = ", num_iterations 
        
        call HYPRE_ParCSRHybridGetNumIterati (gr_hypreSolver, num_iterations, ierr)
!!$        print*, "HYPRE HYBRID (DSC)Num Iterations = ", num_iterations         
        
        call HYPRE_ParCSRHybridGetFinalRelat (gr_hypreSolver, final_res_norm, ierr) 
!!$        print*, "HYPRE HYBRID Relative Residual Norm = ", final_res_norm
        
     end select
     
     if (num_iterations >= gr_hypreMaxIter) then     

        print*,"[gr_hypreSolve]: Nonconvergence in subroutine ", & 
             "gr_hypreSolve after max iterations, final_res_norm =",final_res_norm     
        
     else if (gr_hyprePrintSolveInfo) then
        
        print*, "HYPRE SOLVE: Num Iterations = ", num_iterations 
        print*, "HYPRE SOLVE: Relative Residual Norm = ", final_res_norm          
        
     end if
     
  end if
  

  
  call Timers_stop("gr_hypreSolve") 
  
  return
  
end subroutine gr_hypreSolve
