!!****if* source/physics/materialProperties/MassDiffusivity/MassDiffusivityMain/Constant/MassDiffusivity_init
!!
!! NAME
!!
!!  MassDiffusivity_init
!!
!! SYNOPSIS
!!
!!  MassDiffusivity_init()
!!
!! DESCRIPTION
!!
!!  Initializes MassDiffusivity parameters in MassDiffusivity_data for a 
!!  constant mass diffusivity.
!!
!! ARGUMENTS
!!
!!
!!
!!
!!***

subroutine MassDiffusivity

  use MassDiffusivity_data, ONLY: massdiff_spec_D
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

  

  call RuntimeParameters_get("diff_spec_D", massdiff_spec_D)

end subroutine MassDiffusivity_init
