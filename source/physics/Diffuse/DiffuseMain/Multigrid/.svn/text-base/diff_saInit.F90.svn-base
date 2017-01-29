!!****if* source/physics/Diffuse/DiffuseMain/Multigrid/diff_saInit
!!
!! NAME
!!
!!  diff_saInit
!!
!!
!! SYNOPSIS
!!
!!  call diff_saInit()
!!
!! Description
!!
!!  Initializes local data for Unit Diffuse defined in Module diff_saData.
!!  All the variables here are initialized by calling the
!!  RuntimeParameters_get subroutine. These data variables are for
!!  Unit Scope ->  Diffuse.
!!
!! ARGUMENTS
!!
!!  none  
!!
!! PARAMETERS
!!
!!    diff_scaleFactThermSaTempDiff
!!        factor by which the solution is scaled.
!!    diff_scaleFactThermSaTime
!!        factor by which diffusion time step is scaled.
!!***

subroutine diff_saInit

  use RuntimeParameters_interface, ONLY: RuntimeParameters_get
  use diff_saData, ONLY: diff_scaleFactThermSaTempDiff, diff_scaleFactThermSaTime

  implicit none

  call RuntimeParameters_get('diff_scaleFactThermSaTempDiff',diff_scaleFactThermSaTempDiff)
  call RuntimeParameters_get('diff_scaleFactThermSaTime',diff_scaleFactThermSaTime)

end subroutine diff_saInit
