!!****if* source/Grid/GridStructures/gr_sbInit
!!
!! NAME
!!
!!  gr_sbInit
!!
!! SYNOPSIS
!!
!!  gr_sbInit()
!!
!! DESCRIPTION
!!
!!  Called from Grid_init. Allocate and populate the data structure that holds all information
!!  about the solid bodies. Currently a stub, proper initialization of  the solid bodies 
!!  required in a real application. 
!!
!! ARGUMENTS
!!
!!
!!***


Subroutine gr_sbInit()
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, &
       gr_sbPtNumX, gr_sbPtNumY, gr_sbPtNumZ, gr_sbDebug
!  use Grid_interface, ONLY : Grid_sbBroadcastParticles,Grid_sbSelectMaster
!  use gr_sbInterface, ONLY: gr_sbCreateParticles

  implicit none
  nullify(gr_sbBodyInfo)
  gr_sbNumBodies = 0
  gr_sbPtNumX = 0
  gr_sbPtNumY = 0
  gr_sbPtNumZ = 0
  gr_sbDebug = .true.

!  call Grid_sbSelectMaster()
!  call gr_sbCreateParticles()
!  call Grid_sbBroadcastParticles()


End Subroutine gr_sbInit
