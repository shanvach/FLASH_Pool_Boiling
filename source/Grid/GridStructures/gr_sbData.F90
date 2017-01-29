!!****if* source/Grid/GridStructures/gr_sbData
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
   integer :: comm
   integer :: myPE
   integer :: numProcs
   integer :: bodyMaster
end type solid_body

type(solid_body), save, dimension(:), pointer :: gr_sbBodyInfo
integer, save :: gr_sbNumBodies

!Number of particles to generate within each solid body.
integer, save :: gr_sbPtNumX
integer, save :: gr_sbPtNumY
integer, save :: gr_sbPtNumZ

logical, save :: gr_sbDebug

End Module gr_sbData
