!!****if* source/Grid/GridSolvers/HYPRE_KPDa/gr_hypreData
!!
!! NAME
!!
!!  gr_hypreData
!!
!! SYNOPSIS
!!  use gr_hypreData
!!
!! DESCRIPTION
!!
!!  Defines and stores local data for the HYPRE implementation
!!
!!  
!!***

#include "Flash.h"

Module gr_hypreData 
  
  use gr_interfaceTypeDecl
  
  implicit none

  
  integer, parameter :: range=SELECTED_INT_KIND(16)

  integer,save  :: gr_hypreSolverType
  integer,save  :: gr_hyprePcType 
  real,   save  :: gr_hypreRelTol
  integer,save  :: gr_hypreMaxIter
  integer,save  :: gr_hypreInfoLevel
  logical, save :: gr_hyprePrintSolveInfo
  
  integer (KIND=range), save :: gr_hypreSolver
  integer (KIND=range), save :: gr_hyprePC
  integer (KIND=range), save :: gr_hypreGrid
  integer (KIND=range), save :: gr_hypreStencil
  integer (KIND=range), save :: gr_hypreGraph
  integer (KIND=range), save :: gr_hypreMatA
  integer (KIND=range), save :: gr_hypreVecB
  integer (KIND=range), save :: gr_hypreVecX
  
  integer, allocatable, save ::  gr_hypreLower(:,:)
  integer, allocatable, save ::  gr_hypreUpper(:,:)
  integer, allocatable, save ::  gr_hypreNeghLevels(:,:,:,:) 
  
  integer, save :: gr_hypreNParts
  integer, save :: gr_hypreNVars  
  logical, save :: gr_hypreSetup
  
  integer, save :: gr_hypreRefineMIN
  integer, save :: gr_hypreRefineMAX

  real,    save :: gr_asol, gr_speedlt 

  logical, save :: gr_hypreUseFloor
  real,    save :: gr_hypreFloor
  logical, save :: gr_hypreUse2Norm
  
  type (AllBlockRegions_t), allocatable, save :: gr_hypreSurrBlkSum(:)

  integer, save :: gr_hypreNStep
  
  
end Module gr_hypreData
