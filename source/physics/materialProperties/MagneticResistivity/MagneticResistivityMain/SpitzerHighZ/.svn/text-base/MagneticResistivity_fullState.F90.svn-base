!!****if* source/physics/materialProperties/MagneticResistivity/MagneticResistivityMain/SpitzerHighZ/MagneticResistivity_fullState
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
  use MagneticResistivity_interface, ONLY: MagneticResistivity
  use MagneticResistivity_data, ONLY: mag_useMagneticResistivity, &
       res_mele, res_qele, res_navo, res_speedlt
  use MagneticResistivity_data, ONLY: res_mUnit
  use MagneticResistivity_data, ONLY: res_coef
  use Eos_interface, ONLY: Eos_getAbarZbar

  implicit none

  real, intent(in)  :: solnVec(NUNK_VARS)
  real, intent(out) :: resPar
  real, intent(out) :: resPerp

  real :: dens
  real :: tele
  real :: tion
  real :: nele
  real :: nion
  real :: abar
  real :: zbar
  real :: eqtime

  call Eos_getAbarZbar(solnVec=solnVec,abar=abar,zbar=zbar)

  dens = solnVec(DENS_VAR)
  nion = dens * res_navo / abar
  nele = zbar * nion

  tele = solnVec(TELE_VAR)
  tion = solnVec(TION_VAR)

  call res_ieEquilTime(zbar, abar, tele, tion, nion, eqtime)

  resPerp = res_coef * res_mele / (res_qele**2 * nele * eqtime)

  ! This formula is only valid when the magnetic field is strong (and
  ! only for hydrogen):
  resPar = 1.96 * resPerp

  !! Scale resistivity depending on units
  if (res_mUnit == "SI" .or. res_mUnit == "si" ) then
     resPerp = resPerp*1.e7/(4.*PI)
     resPar  = resPar *1.e7/(4.*PI)
  elseif (res_mUnit == "CGS" .or. res_mUnit == "cgs" ) then
     resPerp = resPerp*res_speedlt**2/(4.*PI)
     resPar  = resPar *res_speedlt**2/(4.*PI)
  else !no unit
     ! Do nothing
     resPerp = resPerp*res_speedlt**2
     resPar  = resPar *res_speedlt**2
  end if
end subroutine MagneticResistivity_fullState
