!!****if* source/physics/Diffuse/DiffuseMain/Multigrid/diff_saData
!!
!! NAME
!!
!!  diff_saData
!!
!! SYNOPSIS
!!  use diff_saData
!!
!! DESCRIPTION
!!
!!  Defines and stores local data for the Multigrid implementation of the main subunit of unit Diffuse.
!!  All the variables defined here should be initialized in Diffuse_init() by calling
!!  RuntimeParameters_get subroutine or similar.
!!  
!!***


Module diff_saData

  use Grid_interface, ONLY: GRID_PDE_BND_DIRICHLET

  implicit none

  logical, save :: updateDiffuse = .TRUE.
  integer, save :: diff_boundary = GRID_PDE_BND_DIRICHLET

  real, save :: diff_scaleFactThermSaTempDiff, diff_scaleFactThermSaTime

end Module diff_saData
