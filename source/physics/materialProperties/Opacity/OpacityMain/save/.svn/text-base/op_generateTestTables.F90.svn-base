!!****if* source/physics/materialProperties/Opacity/OpacityMain/save/op_generateTestTables
!!
!! NAME
!!
!!  op_generateTestTables
!!
!! SYNOPSIS
!!
!!  call op_generateTestTables (character (in) :: tableKind)
!!
!! DESCRIPTION
!!
!!  Generates tabulated Opacities for each species for a specific kind
!!  of table. Currently only the IONMIX kind can be generated.
!!
!!  The names of the files generated will be:
!!
!!                      OPACITIES_xxx.<ext>
!!
!!  where 'xxx' denotes the species number (# of species > 999 induce an error) and
!!  'ext' labels the table kind. Currently only one label is used:
!!
!!                            imx  -->  IONMIX files
!!
!! ARGUMENTS
!!
!!  tableKind     : the kind of tabulated Opacity file where data is going to be read
!!  totalSpecies  : the total number of species.
!!  ngroupsEnergy : the total number of energy groups.
!!
!!***
subroutine op_generateTestTables (tableKind)

  use Driver_interface,     ONLY : Driver_abortFlash
  use Opacity_dataUnitTest, ONLY : op_totalSpecies,          &
                                   op_ngroupsEnergy,         &
                                   op_nstepsTemperature,     &
                                   op_nstepsDensity,         &
                                   op_temperatureFirst,      &
                                   op_temperatureLast,       &
                                   op_temperatureStep,       &
                                   op_ionNumberDensityFirst, &
                                   op_ionNumberDensityLast,  &
                                   op_ionNumberDensityStep,  &
                                   op_log10opacityFirstPA,   &
                                   op_log10opacityFirstPE,   &
                                   op_log10opacityFirstRO,   &
                                   op_log10opacityStepPA,    &
                                   op_log10opacityStepPE,    &
                                   op_log10opacityStepRO

  use op_interface,         ONLY : op_generateTestIonmixTable

  implicit none

  character (len=*), intent (in) :: tableKind

  character (len= 1) :: tabledot  = '.'
  character (len=10) :: tableBase = 'OPACITIES_'
  character (len= 3) :: tableLabel
  character (len=17) :: tableName

  integer :: nstepsDensity
  integer :: nstepsTemperature
  integer :: species

  real    :: densityFirst
  real    :: densityLast
  real    :: densityStep
  real    :: log10opacityFirstPA
  real    :: log10opacityFirstPE
  real    :: log10opacityFirstRO
  real    :: log10opacityStepPA
  real    :: log10opacityStepPE
  real    :: log10opacityStepRO
  real    :: temperatureFirst
  real    :: temperatureLast
  real    :: temperatureStep
!
!
!   ...Branch according to table kind and loop over all species. 
!
!
  if (tableKind == "IONMIX") then

      tableLabel = 'imx'

      do species = 1,op_totalSpecies

         write (tableName,'(A10,I3.3,A1,A3)') tableBase,species,tabledot,tableLabel

         nstepsTemperature   = op_nstepsTemperature     (species)
         nstepsDensity       = op_nstepsDensity         (species)
         temperatureFirst    = op_temperatureFirst      (species)
         temperatureStep     = op_temperatureStep       (species)
         densityFirst        = op_ionNumberDensityFirst (species)
         densityStep         = op_ionNumberDensityStep  (species)
         log10opacityFirstPA = op_log10opacityFirstPA   (species)
         log10opacityFirstPE = op_log10opacityFirstPE   (species)
         log10opacityFirstRO = op_log10opacityFirstRO   (species)
         log10opacityStepPA  = op_log10opacityStepPA    (species)
         log10opacityStepPE  = op_log10opacityStepPE    (species)
         log10opacityStepRO  = op_log10opacityStepRO    (species)

         call op_generateTestIonmixTable (tableName,           &
                                          nstepsTemperature,   &
                                          nstepsDensity,       &
                                          op_ngroupsEnergy,    &
                                          temperatureFirst,    &
                                          temperatureStep,     &
                                          densityFirst,        &
                                          densityStep,         &
                                          log10opacityFirstPA, &
                                          log10opacityFirstPE, &
                                          log10opacityFirstRO, &
                                          log10opacityStepPA,  &
                                          log10opacityStepPE,  &
                                          log10opacityStepRO,  &

                                              temperatureLast, &
                                              densityLast      )

         op_temperatureLast      (species) = temperatureLast
         op_ionNumberDensityLast (species) = densityLast

      end do
  else
      call Driver_abortFlash ('[op_generateTestTables] ERROR: Opacity table kind not recognized')
  end if
!
!
!    ...Ready!
!
!
  return
end subroutine op_generateTestTables
