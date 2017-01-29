!!****if* source/physics/materialProperties/Opacity/OpacityMain/save/Opacity_finalizeUnitTest
!!
!! NAME
!!
!!  Opacity_finalizeUnitTest
!!
!! SYNOPSIS
!!
!!  call Opacity_finalizeUnitTest ()
!!
!! DESCRIPTION
!!
!!  Clean up the unit test opacity runs.
!!
!! ARGUMENTS
!!
!!***
subroutine Opacity_finalizeUnitTest ()

  use Opacity_dataUnitTest

  implicit none
!
!
!    ...Remove the allocated arrays.
!
!
  deallocate (op_nstepsTemperature)
  deallocate (op_nstepsDensity)
  deallocate (op_temperatureFirst)
  deallocate (op_temperatureLast)
  deallocate (op_temperatureStep)
  deallocate (op_massDensityFirst)
  deallocate (op_massDensityLast)
  deallocate (op_massDensityStep)
  deallocate (op_ionNumberDensityFirst)
  deallocate (op_ionNumberDensityLast)
  deallocate (op_ionNumberDensityStep)
  deallocate (op_log10opacityFirstPA)
  deallocate (op_log10opacityFirstPE)
  deallocate (op_log10opacityFirstRO)
  deallocate (op_log10opacityStepPA)
  deallocate (op_log10opacityStepPE)
  deallocate (op_log10opacityStepRO)
  deallocate (op_absorptionKind)
  deallocate (op_emissionKind)
  deallocate (op_transportKind)
  deallocate (op_speciesWeights)
  deallocate (op_massFractions)
!
!
!    ...Ready!
!
!
end subroutine Opacity_finalizeUnitTest
