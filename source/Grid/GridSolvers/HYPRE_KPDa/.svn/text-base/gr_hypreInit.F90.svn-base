!!****if* source/Grid/GridSolvers/HYPRE_KPDa/gr_hypreInit
!!
!! NAME
!!
!!  gr_hypreInit
!!
!!
!! SYNOPSIS
!!
!!  call gr_hypreInit()
!!
!! Description
!!
!!  Initializes local data for Unit HYPRE defined in Module gr_hypreData.
!!  All the variables here are initialized by calling the
!!  RuntimeParameters_get subroutine. These data variables are for
!!  Unit Scope ->  HYPRE.
!!
!! ARGUMENTS
!!
!!  none  
!!***

subroutine gr_hypreInit()
  
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
                                          RuntimeParameters_mapStrToInt
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use gr_hypreData, ONLY : gr_asol, gr_speedlt, gr_hypreSolverType,            &
                           gr_hyprePcType, gr_hypreRelTol, gr_hypreMaxIter,    & 
                           gr_hypreSetup, gr_hyprePrintSolveInfo,              &
                           gr_hypreInfoLevel, gr_hypreUseFloor, gr_hypreFloor, &
                           gr_hypreUse2Norm
  
  implicit none 

#include "constants.h"
 
  character(len=MAX_STRING_LENGTH) :: solver
  character(len=MAX_STRING_LENGTH) :: pc
  
  real :: ssol
  
  call PhysicalConstants_get("speed of light",gr_speedlt)
  call PhysicalConstants_get("Stefan-Boltzmann",ssol)
  
  gr_asol    = 4.0 * ssol / gr_speedlt
  
  gr_hypreSolverType   = HYPRE_PCG
  gr_hyprePcType       = HYPRE_AMG
  gr_hypreRelTol       = 1.0E-12
  gr_hypreMaxIter      = 10000
  
  call RuntimeParameters_get("gr_hypreSolverType", solver)
  call RuntimeParameters_get("gr_hyprePcType", pc)
  
  call RuntimeParameters_mapStrToInt(solver,gr_hypreSolverType)
  call RuntimeParameters_mapStrToInt(pc,gr_hyprePcType)

  call RuntimeParameters_get("gr_hypreRelTol",  gr_hypreRelTol)
  call RuntimeParameters_get("gr_hypreMaxIter", gr_hypreMaxIter)

  call RuntimeParameters_get("gr_hyprePrintSolveInfo", gr_hyprePrintSolveInfo)
  
  call RuntimeParameters_get("gr_hypreInfoLevel", gr_hypreInfoLevel)
  
  call RuntimeParameters_get("gr_hypreUseFloor", gr_hypreUseFloor)
  call RuntimeParameters_get("gr_hypreFloor",    gr_hypreFloor)

  call RuntimeParameters_get("gr_hypreUse2Norm",    gr_hypreUse2Norm)
  
  gr_hypreSetup = .FALSE.   
  
  !! setup the HYPRE diffusion solver.
  call gr_hypreSetupSolver ()
  
  
end subroutine gr_hypreInit
