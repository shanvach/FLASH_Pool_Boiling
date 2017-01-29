!!****if* source/physics/materialProperties/Opacity/localAPI/save/op_generateTestIonmixTable
!!
!! NAME
!!
!!  op_generateTestIonmixTable
!!
!! SYNOPSIS
!!
!!  call op_generateTestIonmixTable (character (in)  :: tableName,
!!                                   integer   (in)  :: nstepsTemperature,
!!                                   integer   (in)  :: nstepsDensity,
!!                                   integer   (in)  :: ngroupsEnergy,
!!                                   real      (in)  :: temperatureFirst,
!!                                   real      (in)  :: temperatureStep,
!!                                   real      (in)  :: densityFirst,
!!                                   real      (in)  :: densityStep,
!!                                   real      (in)  :: log10opacityFirstPA,
!!                                   real      (in)  :: log10opacityFirstPE,
!!                                   real      (in)  :: log10opacityFirstRO,
!!                                   real      (in)  :: log10opacityStepPA,
!!                                   real      (in)  :: log10opacityStepPE,
!!                                   real      (in)  :: log10opacityStepRO,
!!                                   real      (out) :: temperatureLast,
!!                                   real      (out) :: densityLast)
!!
!! DESCRIPTION
!!
!!  This routine creates an opacities datafile in the IONMIX output format. This
!!  routine is meant for checking correctness of the Opacity unit and should never
!!  be used in real applications. The following picture is a layout of the table
!!  for one energy group g:
!!
!!
!!                             Df --> Dstep
!!                             ----------------------------
!!                         Tf |Of --> Ostep                |
!!                            |                            |
!!                         |  | |                          |
!!                         |  | |                          |
!!                         v  | v                          |
!!                            |                            |
!!                      Tstep | Ostep                      |
!!                            |                            |
!!                            |                            |
!!                            |                            |
!!                            |                            |
!!                            |                            |
!!                            |                            |
!!                            |                            |
!!                            |                            |
!!                            |                            |
!!                            |                            |
!!                            |                          Ol|
!!                             ----------------------------
!!
!!
!!  The values of the three first-step pairs for the temperature, density and opacities
!!  are supplied in argument.
!!
!!  The opacity values generated for the three kind of opacities are then defined
!!  as follows on the temperature/density/energy grid:
!!
!!            opacity values -->  O (t,d,g) = 10 ** {g * [log10 (Of) + (t + d - 2) * log10 (Ostep)]}
!!
!!  where 't' and 'd' denote coordinate counts on the  temperature and density grid.
!!
!!  Note the use of the powers of 10. This is done to obtain the intended opacity tables when
!!  using the 'useLogTables' = .true. option. The energy groups are simply given the (real)
!!  values of their indices g = 1,ngroupsEnergy
!!
!! ARGUMENTS
!!
!!  tableName           : the name of the created IONMIX file
!!  nstepsTemperature   : defines the temperature grid
!!  nstepsDensity       : defines the density grid
!!  ngroupsEnergy       : defines the energy grid
!!  temperatureFirst    : the initial temperature on the temperature grid
!!  temperatureStep     : the   step  temperature on the temperature grid
!!  densityFirst        : the initial density on the density grid
!!  densityStep         : the   step  density on the density grid
!!  log10opacityFirstPA : the initial log10 Planck Absorption opacity on the (temperature,density) grid
!!  log10opacityFirstPE : the initial log10 Planck  Emission  opacity on the (temperature,density) grid
!!  log10opacityFirstRO : the initial log10      Rosseland    opacity on the (temperature,density) grid
!!  log10opacityStepPA  : the   step  log10 Planck Absorption opacity on the (temperature,density) grid
!!  log10opacityStepPE  : the   step  log10 Planck  Emission  opacity on the (temperature,density) grid
!!  log10opacityStepRO  : the   step  log10      Rosseland    opacity on the (temperature,density) grid
!!  temperatureLast     : the   last  temperature on the temperature grid
!!  densityLast         : the   last  density on the density grid
!!
!!***
subroutine op_generateTestIonmixTable (tableName,            &
                                       nstepsTemperature,    &
                                       nstepsDensity,        &
                                       ngroupsEnergy,        &
                                       temperatureFirst,     &
                                       temperatureStep,      &
                                       densityFirst,         &
                                       densityStep,          &
                                       log10opacityFirstPA,  &
                                       log10opacityFirstPE,  &
                                       log10opacityFirstRO,  &
                                       log10opacityStepPA,   &
                                       log10opacityStepPE,   &
                                       log10opacityStepRO,   &
                                            temperatureLast, &
                                            densityLast      )

  implicit none

  character (len=*), intent (in)  :: tableName
  integer,           intent (in)  :: ngroupsEnergy
  integer,           intent (in)  :: nstepsDensity
  integer,           intent (in)  :: nstepsTemperature
  real,              intent (in)  :: densityFirst
  real,              intent (out) :: densityLast
  real,              intent (in)  :: densityStep
  real,              intent (in)  :: log10opacityFirstPA
  real,              intent (in)  :: log10opacityFirstPE
  real,              intent (in)  :: log10opacityFirstRO
  real,              intent (in)  :: log10opacityStepPA
  real,              intent (in)  :: log10opacityStepPE
  real,              intent (in)  :: log10opacityStepRO
  real,              intent (in)  :: temperatureFirst
  real,              intent (out) :: temperatureLast
  real,              intent (in)  :: temperatureStep

  return
end subroutine op_generateTestIonmixTable
