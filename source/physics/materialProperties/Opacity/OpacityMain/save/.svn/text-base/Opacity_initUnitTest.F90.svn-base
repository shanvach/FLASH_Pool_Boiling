!!****if* source/physics/materialProperties/Opacity/OpacityMain/save/Opacity_initUnitTest
!!
!! NAME
!!
!!  Opacity_initUnitTest
!!
!! SYNOPSIS
!!
!!  call  Opacity_initUnitTest (integer (in):: myPE)
!!
!! DESCRIPTION
!!
!!  This routine initializes the unitTest for testing the Opacity unit. Normally called
!!  from Driver_initFlash within a Simulation unitTest. See for example
!!  the directory 'source/Simulation/SimulationMain/unitTest/Opacity'.
!!
!!  Prepares test opacity input and tables to check on the routines that interpolate
!!  the opacity tables.
!!
!! ARGUMENTS
!!
!!  myPE -- local processor number 
!!
!!***
subroutine Opacity_initUnitTest (myPE)

  use Opacity_dataUnitTest

  use Multispecies_interface,      ONLY : Multispecies_getProperty
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use op_interface,                ONLY : op_generateTestOpacityInput, &
                                          op_generateTestTables

  implicit none

# include "Flash.h"  
# include "Multispecies.h"

  integer, intent (in) :: myPE

  character (len=6) :: tableKind = 'IONMIX'

  integer :: n
  integer :: species

  real    :: a,b
  real    :: speciesWeight
  real    :: ionDensityFirst
  real    :: ionDensityStep
  real    :: massDensityFirst
!
!
!   ...Check the # of species and energy groups. 
!
!
  if (NSPECIES > 0) then
      op_totalSpecies = NSPECIES
  else
      call Driver_abortFlash ('[Opacity_initUnitTest] ERROR: no species found (NSPECIES =< 0)')
  end if

  if (MGD_NGROUPS > 0) then
      op_ngroupsEnergy = MGD_NGROUPS
  else
      call Driver_abortFlash ('[Opacity_initUnitTest] ERROR: no energy groups found (MGD_NGROUPS =< 0)')
  end if
!
!
!    ...Get the external data.
!
!
  call RuntimeParameters_get ("opacity_useLogTables",   op_useLogTables )
  call PhysicalConstants_get ("Avogadro",               op_Avogadro     )
!
!
!   ...Allocate the needed arrays. 
!
!
  allocate (op_speciesWeights        (1:op_totalSpecies))
  allocate (op_massFractions         (1:op_totalSpecies))
  allocate (op_nstepsTemperature     (1:op_totalSpecies))
  allocate (op_nstepsDensity         (1:op_totalSpecies))
  allocate (op_temperatureFirst      (1:op_totalSpecies))
  allocate (op_temperatureLast       (1:op_totalSpecies))
  allocate (op_temperatureStep       (1:op_totalSpecies))
  allocate (op_massDensityFirst      (1:op_totalSpecies))
  allocate (op_massDensityLast       (1:op_totalSpecies))
  allocate (op_massDensityStep       (1:op_totalSpecies))
  allocate (op_ionNumberDensityFirst (1:op_totalSpecies))
  allocate (op_ionNumberDensityLast  (1:op_totalSpecies))
  allocate (op_ionNumberDensityStep  (1:op_totalSpecies))
  allocate (op_log10opacityFirstPA   (1:op_totalSpecies))
  allocate (op_log10opacityFirstPE   (1:op_totalSpecies))
  allocate (op_log10opacityFirstRO   (1:op_totalSpecies))
  allocate (op_log10opacityStepPA    (1:op_totalSpecies))
  allocate (op_log10opacityStepPE    (1:op_totalSpecies))
  allocate (op_log10opacityStepRO    (1:op_totalSpecies))
  allocate (op_absorptionKind        (1:op_totalSpecies))
  allocate (op_emissionKind          (1:op_totalSpecies))
  allocate (op_transportKind         (1:op_totalSpecies))
!
!
!   ...Set the parameters for testing the opacity tables. The goal is to obtain
!      a reasonable opacity table. The temperature scale is in eV and presents
!      no problem. The density introduced first is the mass density for each
!      species from which later on the ion # densities will be calculated.
!      The density scale for the opacity tables should be in ion number densities,
!      i.e. # of ions / cm^3. The ion number densities are huge numbers
!      (~ Avogadro's number) and are dependent on the mass fractions and species
!      weights of each species. Their evaluation must be deferred to at a
!      later stage (see below).
!
!
  do species = 1,op_totalSpecies

     call Multispecies_getProperty (SPECIES_BEGIN - 1 + species , A , speciesWeight)
     write (*,*) ' A in initUnitTest = ',speciesWeight

     op_speciesWeights       (species) = speciesWeight
     op_nstepsTemperature   (species) = 10
     op_nstepsDensity       (species) = 10
     op_temperatureFirst    (species) = ten        ! in eV                              , will be converted to log10 in tables
     op_temperatureStep     (species) = ten        ! takes x10 the previous temperature , will be converted to log10 in tables
     op_massdensityFirst    (species) = one        ! in g / cm^3                        , will be used to form ion # densities
     op_massdensityStep     (species) = two        ! takes x2 the previous mass density , will be used to form ion # densities
     op_log10opacityFirstPA (species) = one        ! in cm^2 / g                        , will be converted to 10 ** in tables
     op_log10opacityFirstPE (species) = one        ! in cm^2 / g                        , will be converted to 10 ** in tables
     op_log10opacityFirstRO (species) = one        ! in cm^2 / g                        , will be converted to 10 ** in tables
     op_log10opacityStepPA  (species) = one        ! takes x10 the previous opacity     , will be converted to 10 ** in tables
     op_log10opacityStepPE  (species) = 2 * one    ! takes x100 the previous opacity    , will be converted to 10 ** in tables
     op_log10opacityStepRO  (species) = 3 * one    ! takes x1000 the previous opacity   , will be converted to 10 ** in tables

     n = mod (species,3)

     if (n == 1) then

         op_absorptionKind (species) = "Planck Absorption"
         op_emissionKind   (species) = "Planck Absorption"
         op_transportKind  (species) = "Planck Absorption"

     else if (n == 2) then

         op_absorptionKind (species) = "Planck Emission"
         op_emissionKind   (species) = "Planck Emission"
         op_transportKind  (species) = "Planck Emission"

     else

         op_absorptionKind (species) = "Rosseland"
         op_emissionKind   (species) = "Rosseland"
         op_transportKind  (species) = "Rosseland"

     end if

  end do
!
!
!   ...Set the mass fractions for each species. This is done such that the species
!      with the largest species weight occurs least.
!
!
  a = one / real (op_totalSpecies - 1)
  b = one / sum (op_speciesWeights (:))

  op_massFractions (:) = a * (one - b * op_speciesWeights (:))
!
!
!   ...We are now ready to evaluate the ion # densities for each species.
!
!
  op_ionNumberDensityFirst (:) = op_Avogadro * op_massFractions (:) * op_massDensityFirst (:) / op_speciesWeights (:)
  op_ionNumberDensityStep  (:) = op_massDensityStep (:)
!
!
!   ...Generate the tables and the required input. 
!
!
  call op_generateTestOpacityInput ()
  call op_generateTestTables       (tableKind)
!
!
!    ...Ready!
!
!
  return
end subroutine Opacity_initUnitTest
