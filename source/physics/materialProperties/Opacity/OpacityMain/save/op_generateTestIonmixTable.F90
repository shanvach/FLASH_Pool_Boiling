!!****if* source/physics/materialProperties/Opacity/OpacityMain/save/op_generateTestIonmixTable
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

  use Driver_interface,  ONLY : Driver_abortFlash

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

  character (len=50) :: dummyLine = ' This is a mockup IONMIX file for testing purposes'

  logical :: fileExists

  integer :: fileUnit
  integer :: n,t,d,g
  integer :: notneededData

  real    :: dummyData = 0.0
  real    :: log10Density
  real    :: log10DensityStep
  real    :: log10DensityFirst
  real    :: log10Temperature
  real    :: log10TemperatureStep
  real    :: log10TemperatureFirst
!
!
!   ...Check the supplied temperature, density and energy grids.
!
!
  if (nstepsTemperature <= 0) then
      call Driver_abortFlash ('[op_generateTestIonmixTable] ERROR: Bad temperature grid')
  end if

  if (nstepsDensity <= 0) then
      call Driver_abortFlash ('[op_generateTestIonmixTable] ERROR: Bad density grid')
  end if

  if (ngroupsEnergy <= 0) then
      call Driver_abortFlash ('[op_generateTestIonmixTable] ERROR: Bad energy grid')
  end if
!
!
!   ...Check and open the IONMIX file.
!
!
  open (unit = fileUnit,   &
        file = tableName,  &
        form = 'formatted' )
!
!
!   ...Write the temperature, density and energy group grids in IONMIX style.
!
!
  write (fileUnit,'(2I10)') nstepsTemperature , nstepsDensity
  write (fileUnit,'(A50)')  dummyLine
  write (fileUnit,'(A50)')  dummyLine

  log10DensityStep      = log10 (densityStep)
  log10DensityFirst     = log10 (densityFirst)
  log10TemperatureStep  = log10 (temperatureStep)
  log10TemperatureFirst = log10 (temperatureFirst)

  write (fileUnit,'(4E12.6,I12)') log10DensityStep,       &
                                  log10DensityFirst,      &
                                  log10TemperatureStep,   &
                                  log10TemperatureFirst,  &
                                  ngroupsEnergy
!
!
!   ...Calculate the maximum temperature and density on the grid. These values
!      will be returned via argument and stored for further reference.
!
!
  temperatureLast = 10 ** (log10TemperatureFirst + real (nstepsTemperature - 1) * log10TemperatureStep)
  densityLast     = 10 ** (log10DensityFirst     + real (nstepsDensity     - 1) * log10DensityStep)
!
!
!   ...Write placeholders (zeros) for not needed data in the IONMIX file.
!
!
  notneededData =   nstepsDensity * nstepsTemperature

  write (fileUnit,'(4E12.6)') (dummyData, n = 1,notneededData)
  write (fileUnit,'(4E12.6)') (dummyData, n = 1,notneededData)
  write (fileUnit,'(4E12.6)') (dummyData, n = 1,notneededData)
  write (fileUnit,'(4E12.6)') (dummyData, n = 1,notneededData)
!
!
!   ...Write the energy groups.
!
!
  write (fileUnit,'(4E12.6)') (real (g), g = 1,ngroupsEnergy+1)
!
!
!   ...Write the Rosseland, Planck Absorption and Planck Emission opacities in that
!      order (same order as in  IONMIX tables).
!
!
  write (fileUnit,'(4E12.6)') (((   &

         10 ** ( real (g) * (log10opacityFirstRO + real (t+d-2) * log10opacityStepRO) ) &

       , t = 1,nstepsTemperature) &
       , d = 1,nstepsDensity    ) &
       , g = 1,ngroupsEnergy    )

  write (fileUnit,'(4E12.6)') (((   &

         10 ** ( real (g) * (log10opacityFirstPA + real (t+d-2) * log10opacityStepPA) ) &

       , t = 1,nstepsTemperature) &
       , d = 1,nstepsDensity    ) &
       , g = 1,ngroupsEnergy    )


  write (fileUnit,'(4E12.6)') (((   &

         10 ** ( real (g) * (log10opacityFirstPE + real (t+d-2) * log10opacityStepPE) ) &

       , t = 1,nstepsTemperature) &
       , d = 1,nstepsDensity    ) &
       , g = 1,ngroupsEnergy    )
!
!
!   ...Close the generated IONMIX file.
!
!
  close (fileUnit)
!
!
!   ...Ready! 
!
!
  return
end subroutine op_generateTestIonmixTable
