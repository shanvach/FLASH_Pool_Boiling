!!****if* source/Grid/GridStructures/withTriangles/gr_sbData
!!
!! NAME
!!  gr_sbData
!!
!! SYNOPSIS
!!
!!  use gr_sbData
!!
!!
!!***

#include "constants.h"

Module gr_sbData
implicit none

type solid_body
   real, pointer, dimension(:,:) :: particles
   real, dimension(2,MDIM) :: boundBox
   real, pointer, dimension(:,:) :: triangleCentroids
   integer :: comm
   integer :: myPE
   integer :: numProcs
   integer :: bodyMaster
   real, allocatable, dimension(:) :: xb, yb, zb !Vertex coordinates
end type solid_body

real, allocatable, dimension(:,:) :: aelem !Triangle elements

type(solid_body), save, dimension(:), pointer :: gr_sbBodyInfo
integer, save :: gr_sbNumBodies

integer, save :: NumAelem
integer, save :: NumVertices

!Number of particles to generate within each solid body.
integer, save :: gr_sbPtNumX
integer, save :: gr_sbPtNumY
integer, save :: gr_sbPtNumZ

integer, save :: totalPart

logical, save :: gr_sbDebug

End Module gr_sbData
