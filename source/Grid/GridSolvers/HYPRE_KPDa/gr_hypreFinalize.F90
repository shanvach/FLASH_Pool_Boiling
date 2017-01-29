!!****if* source/Grid/GridSolvers/HYPRE_KPDa/gr_hypreFinalize
!!
!! NAME
!!
!!  gr_hypreFinalize
!!
!!
!! SYNOPSIS
!!
!!  call gr_hypreFinalize ()
!!
!! Description
!!  
!!  Release HYPRE objects and other temporary storage from memory.
!!
!! ARGUMENTS
!!
!!  none  
!!
!! PARAMETERS
!!
!!***

subroutine gr_hypreFinalize()
  
  implicit none
  
  !! Destroy the HYPRE solver object.
  call gr_hypreDestroySolver()
  call gr_hypreDestroyGrid()
  
end subroutine gr_hypreFinalize
