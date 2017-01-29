!!****if* source/physics/materialProperties/MassDiffusivity/MassDiffusivityMain/Constant/MassDiffusivity
!!
!! NAME
!!
!!  MassDiffusivity
!!
!! SYNOPSIS
!!
!!  MassDiffusivity(real,INTENT(in)  :: xtemp,
!!                  real,INTENT(in)  :: xden,
!!                  real,INTENT(in)  :: massfrac,
!!                  real,INTENT(out)  :: diffusivity)
!!
!! DESCRIPTION
!!
!!  A generic mass diffusivity routine.  Returns a constant diffusivity from 
!!  the runtime parameter, 'diff_spec_D'.
!!
!! ARGUMENTS
!!
!!   xtemp :  temperature in K 
!!
!!   xden : density in g/cm**3
!!
!!   massfrac : mass fractions of the composition
!!
!!   diffusivity : mass diffusivity
!!
!!
!!
!!***


subroutine MassDiffusivity(xtemp,xden,massfrac,diffusivity)
  
  use MassDiffusivity_data, ONLY : massdiff_spec_D

  implicit none
#include "constants.h"
#include "Flash.h"

  real,INTENT(in)    :: xtemp
  real,INTENT(in)    :: xden
  real,INTENT(in)    :: massfrac(NSPECIES)
  real,INTENT(out) :: diffusivity

  diffusivity = massdiff_spec_D

  return 
end subroutine MassDiffusivity


