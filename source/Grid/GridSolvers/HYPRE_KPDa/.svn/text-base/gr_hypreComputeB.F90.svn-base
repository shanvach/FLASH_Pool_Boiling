!!****if* source/Grid/GridSolvers/HYPRE_KPDa/gr_hypreComputeB
!!
!!  NAME 
!!
!!  gr_hypreComputeB
!!
!!  SYNOPSIS
!!
!!  call gr_hypreComputeB (integer, intent(IN) :: iVar
!!                         integer, intent(IN) :: iFactorB
!!                         integer, intent(IN) :: iFactorA
!!                         integer, OPTIONAL, intent(IN) :: iFactorD
!!                         real, intent(IN) :: dt
!!                         real, intent(IN) :: theta
!!                         integer, intent(IN) :: blockCount
!!                         integer,dimension(blockCount),intent(IN) :: blockList
!!                         integer, intent(IN) :: bcTypes(6)
!!                         real,    intent(IN) :: bcValues(2,6))
!!
!!
!!  DESCRIPTION 
!!   Computes the RHS of AX=B using a precomputed matrix M (stored in diff_A) such that B=MX.
!!   A MatVec product is performed to compute B and additional source terms are added 
!!   if required. A lot of the arguments are not required by this routine but are nevertheless
!!   provided so that B can be computed outside of HYPRE (if needed).
!!  
!!
!! ARGUMENTS
!!   iVar          : Variable on which the diffusion operatorion is performed (e.g TEMP_VAR)
!!   blockCount    : The number of blocks in the list.   
!!   blockList     : The list of blocks on which the solution must be updated.   
!!   iFactorA      :| Are factors in the equation with spatial variation. Factor D is  optional 
!!   iFactorB      :| and is generally used to represent emission in MGD. 
!!   iFactorD      :| 
!!   theta         : varies scheme (0-> Explicit, 1-> backward euler, 0.5 -> Crank Nicholson
!!   dt            : The time step (not used).
!!   bcTypes       : Presently OUTFLOW, VACUUM is supported, DIRICHLET is untested (not used).
!!   bcValues      : Values of iVar,iFactorB on boundary (DIRICHLET), (not used).
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES:
!!
!!   Uses HYPRE library.
!!  
!!
!!***

subroutine gr_hypreComputeB (blockCount, blockList, iVar, iFactorA, iFactorB, &
                             dt, theta, bcTypes, bcValues, iFactorD)
  
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_data,        ONLY : gr_geometry
  use Grid_interface,   ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr,    &
                               Grid_getBlkIndexLimits, &
                               Grid_getBlkData
  use gr_hypreData,     ONLY : gr_hypreMatA, gr_hypreVecB, gr_hypreVecX, &
                               gr_hypreLower, gr_hypreRefineMIN, gr_hypreUpper

  implicit none
  
#include "Flash.h"
#include "constants.h" 
  
  integer, intent(IN) :: iVar
  integer, intent(IN) :: iFactorB
  integer, intent(IN) :: iFactorA
  integer, OPTIONAL, intent(IN) :: iFactorD
  real, intent(IN) :: dt
  real, intent(IN) :: theta
  integer, intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList
  integer, intent(IN) :: bcTypes(6)
  real,    intent(IN) :: bcValues(2,6)
  
  integer::  ierr

  call Timers_start("gr_hypreComputeB")
  
  !This is a collective call finalizing the vector assembly
  call HYPRE_SStructVectorAssemble(gr_hypreVecB, ierr)

  if (theta /= 0.0) then !! Implicit
     call HYPRE_SStructMatrixMatvec (-(1.0-theta)/theta, gr_hypreMatA, gr_hypreVecX, 1.0, gr_hypreVecB, ierr)
  else
     call HYPRE_SStructMatrixMatvec (-1, gr_hypreMatA, gr_hypreVecX, 0.0, gr_hypreVecB, ierr)
  end if

  call HYPRE_SStructVectorGather(gr_hypreVecB, ierr)    
  

  
  call Timers_stop("gr_hypreComputeB")
  
end subroutine gr_hypreComputeB
