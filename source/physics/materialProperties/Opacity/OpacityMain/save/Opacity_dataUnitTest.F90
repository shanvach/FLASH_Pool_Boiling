!!****if* source/physics/materialProperties/Opacity/OpacityMain/save/Opacity_dataUnitTest
!!
!! NAME
!!
!!  Opacity_dataUnitTest
!!
!! SYNOPSIS
!!  use Opacity_dataUnitTest
!!
!! DESCRIPTION
!!
!!  Defines and stores the local data for the opacity unit test.
!!  
!!***
module Opacity_dataUnitTest
  
  implicit none

  logical, save :: op_useLogTables

  integer, save :: op_ngroupsEnergy
  integer, save :: op_totalSpecies

  real,    save :: op_Avogadro

  real, parameter :: zero =  0.0
  real, parameter :: one  =  1.0
  real, parameter :: two  =  2.0
  real, parameter :: ten  = 10.0

  character (len=20), allocatable, save :: op_absorptionKind        (:)
  character (len=20), allocatable, save :: op_emissionKind          (:)
  character (len=20), allocatable, save :: op_transportKind         (:)

  integer,            allocatable, save :: op_nstepsTemperature     (:)
  integer,            allocatable, save :: op_nstepsDensity         (:)
  real,               allocatable, save :: op_temperatureFirst      (:)
  real,               allocatable, save :: op_temperatureLast       (:)
  real,               allocatable, save :: op_temperatureStep       (:)
  real,               allocatable, save :: op_massDensityFirst      (:)
  real,               allocatable, save :: op_massDensityLast       (:)
  real,               allocatable, save :: op_massDensityStep       (:)
  real,               allocatable, save :: op_ionNumberDensityFirst (:)
  real,               allocatable, save :: op_ionNumberDensityLast  (:)
  real,               allocatable, save :: op_ionNumberDensityStep  (:)
  real,               allocatable, save :: op_log10opacityFirstPA   (:)
  real,               allocatable, save :: op_log10opacityFirstPE   (:)
  real,               allocatable, save :: op_log10opacityFirstRO   (:)
  real,               allocatable, save :: op_log10opacityStepPA    (:)
  real,               allocatable, save :: op_log10opacityStepPE    (:)
  real,               allocatable, save :: op_log10opacityStepRO    (:)
  real,               allocatable, save :: op_speciesWeights        (:)
  real,               allocatable, save :: op_massFractions         (:)

end module Opacity_dataUnitTest
