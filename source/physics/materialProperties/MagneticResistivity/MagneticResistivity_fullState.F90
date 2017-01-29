!!****f* source/physics/materialProperties/MagneticResistivity/MagneticResistivity_fullState
!!
!! NAME
!!  MagneticResistivity_fullState
!!
!! SYNOPSIS
!!  call MagneticResistivity_fullState(real(in)    :: solnVec(NUNK_VARS),
!!                            OPTIONAL,real(out)   :: isochoricRes,
!!                            OPTIONAL,real(out)   :: diffCoeff,
!!                            OPTIONAL,integer(in) :: component)
!!
!! DESCRIPTION
!!
!! Computes the Spitzer electron Magnetic Resistivity for all materials,
!! including those with Z > 1. The specific equations used here all
!! come from "The Physics of Inertial Fusion" by Atzeni.
!!
!!  Returns thermal Magnetic Resistivity and/or diffusivity coefficients.
!!
!! ARGUMENTS
!!
!!   solnVec  :   solution state, a vector from UNK with all variables
!!   isochoricRes  :   isochoric Magnetic Resistivity
!!   diffCoeff :   diffusion coefficient ( = isochoricRes/(rho*cv))
!!   component  :   In 3T applications, select component for which MagneticResistivity
!!                  and diffusivity are requested, 1 for ions, 2 for electrons, 
!!                  3 for radiation.
!!
!!***

#include "Flash.h"
#include "constants.h"  

subroutine MagneticResistivity_fullState(solnVec,resPar, resPerp)
  implicit none

  real, intent(in)  :: solnVec(NUNK_VARS)
  real, intent(out) :: resPar
  real, intent(out) :: resPerp

  ! Stub
  resPar  = 0.0
  resPerp = 0.0
end subroutine MagneticResistivity_fullState
