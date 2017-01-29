!!****if* source/physics/Diffuse/DiffuseMain/UG/diff_saFinalize
!!
!! NAME
!!
!!  diff_saFinalize
!!
!!
!! SYNOPSIS
!!
!!  call diff_saFinalize
!!
!! Description
!!
!!
!!
!! ARGUMENTS
!!
!!  none  
!!
!! PARAMETERS
!!
!!***

subroutine diff_saFinalize  
  
  implicit none
  
  !! Destroy the HYPRE solver object.
  call diff_destroyHypreSolver() 
  call diff_destroyHypreGrid()
  
end subroutine diff_saFinalize
