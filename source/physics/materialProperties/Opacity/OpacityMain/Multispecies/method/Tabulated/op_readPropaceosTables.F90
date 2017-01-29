!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Tabulated/op_readPropaceosTables
!!
!! NAME
!!
!!  op_readPropaceosTables
!!
!! SYNOPSIS
!!
!!  call op_readPropaceosTables (character (in) :: tableName (len=80),
!!                              logical   (in) :: needPATable,
!!                              logical   (in) :: needPETable,
!!                              logical   (in) :: needROTable,
!!                              integer   (in) :: indexPA,
!!                              integer   (in) :: indexPE,
!!                              integer   (in) :: indexRO)
!!
!! DESCRIPTION
!!
!!  Reads tabulated opacities from an PROPACEOS datafile output. The tabulated opacities
!!  will be stored into the 4-dimensional arrays:
!!
!!             op_PlanckAbsorptionTables (t,d,g,indexPA)  (in cm^2/g)
!!             op_PlanckEmissionTables   (t,d,g,indexPE)  (in cm^2/g)
!!             op_RosselandTables        (t,d,g,indexRO)  (in cm^2/g)
!!
!!    where:   t = temperature (K) index
!!             d = ion number density (# ions/cm^3) index
!!             g = energy group (eV) index
!!       indexXX = table counting index for Opacity kind XX (XX = PA,PE,RO)
!!
!!  The number of temperature and density indices will be stored in the 1-dimensional arrays: 
!!
!!             op_nstepsDensityPA  (indexPA)
!!             op_nstepsDensityPE  (indexPE)
!!             op_nstepsDensityRO  (indexRO)
!!
!!             op_nstepsTemperaturePA  (indexPA)
!!             op_nstepsTemperaturePE  (indexPE)
!!             op_nstepsTemperatureRO  (indexRO)
!!
!!  The actual temperatures and densities will be stored in the 2-dimensional arrays:
!!
!!           op_tableDensityPA     (d,indexPA) ; d = 1,op_nstepsDensityPA (indexPA)       (in # ions/cm^3)
!!           op_tableDensityPE     (d,indexPE) ; d = 1,op_nstepsDensityPE (indexPE)       (in # ions/cm^3)
!!           op_tableDensityRO     (d,indexRO) ; d = 1,op_nstepsDensityRO (indexRO)       (in # ions/cm^3)
!!
!!           op_tableTemperaturePA (t,indexPA) ; t = 1,op_nstepsTemperaturePA (indexPA)   (in K)
!!           op_tableTemperaturePE (t,indexPE) ; t = 1,op_nstepsTemperaturePA (indexPE)   (in K)
!!           op_tableTemperatureRO (t,indexRO) ; t = 1,op_nstepsTemperaturePA (indexRO)   (in K)
!!
!!  The energy group boundaries will be stored temporarily in the 1-dimensional array:
!!
!!           op_tabulatedEnergyBoundaries (g) ; g = 1,op_nEnergyGroups+1     (in eV)
!!
!!  where
!!           op_tabulatedEnergyBoundaries (g)   = lower boundary of group 'g'
!!           op_tabulatedEnergyBoundaries (g+1) = upper boundary of group 'g'
!!
!!  and will be checked against the energy group boundaries stored for the opacity unit.
!!  A discrepancy within the predefined energy difference tolerance will signal an inconsistency
!!  in the tables generated for the current problem.
!!
!! ARGUMENTS
!!
!!  tableName   : the name of the PROPACEOS file
!!  needPATable : if yes, Planck Absorption Opacities are needed from the PROPACEOS table
!!  needPETable : if yes, Planck   Emission Opacities are needed from the PROPASEOS table
!!  needROTable : if yes,         Rosseland Opacities are needed from the PROPACEOS table
!!  indexPA     : table counting index where Planck Absorption Opacities will be placed
!!  indexPE     : table counting index where Planck   Emission Opacities will be placed
!!  indexRO     : table counting index where Planck  Transport Opacities will be placed
!!
!! NOTES
!!
!!  Since the PROPACEOS tables are produced with temperature units in eV, we must convert these
!!  here into units of K. If the use of logarithmic Tables has been specified, we need to
!!  be careful to substitute the exact zeros of the tables by the smallest representable
!!  positive number in order to avoid NaN's when taking the logarithms.
!!
!!***
subroutine op_readPropaceosTables (tableName,   &
                                 needPATable, &
                                 needPETable, &
                                 needROTable, &
                                 indexPA,     &
                                 indexPE,     &
                                 indexRO      )

  use Driver_interface,  ONLY : Driver_abortFlash

  use Opacity_data,      ONLY : op_nEnergyGroups,             &
                                op_energyGroupBoundaries

  use op_tabulatedData,  ONLY : op_useLogTables,              &
                                op_tableEnergyTolerance,      &
                                op_nstepsDensityPA,           &
                                op_nstepsDensityPE,           &
                                op_nstepsDensityRO,           &
                                op_nstepsTemperaturePA,       &
                                op_nstepsTemperaturePE,       &
                                op_nstepsTemperatureRO,       &
                                op_tabulatedEnergyBoundaries, &
                                op_tableDensityPA,            &
                                op_tableDensityPE,            &
                                op_tableDensityRO,            &
                                op_tableTemperaturePA,        &
                                op_tableTemperaturePE,        &
                                op_tableTemperatureRO,        &
                                op_RosselandTables,           &
                                op_PlanckAbsorptionTables,    &
                                op_PlanckEmissionTables

  use op_numericsData,   ONLY : zero,ten,                     &
                                op_eV2Kelvin,                 &
                                op_smallestPositiveNumber

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Opacity.h"

  character (len=80), intent (in) :: tableName
  logical,            intent (in) :: needPATable
  logical,            intent (in) :: needPETable
  logical,            intent (in) :: needROTable
  integer,            intent (in) :: indexPA
  integer,            intent (in) :: indexPE
  integer,            intent (in) :: indexRO

  character (len=80) :: dummyLine

  logical :: fileExists

  integer :: fileUnit
  integer :: ut_getFreeFileUnit
  integer :: notneededData
  integer :: step
  integer :: n,t,d,g
  integer :: nEnergyGroups
  integer :: nstepsDensity
  integer :: nstepsTemperature

  real    :: dummyData
  real    :: energyDifference
  real    :: newZero

  real, allocatable :: densities(:)
  real, allocatable :: temperatures(:)
!
!
!   ...Check and open the opacity file.
!
!
  inquire (file = tableName , exist = fileExists)

  if (.not.fileExists) then
       call Driver_abortFlash ('[op_readPropaceosTables] ERROR: no PROPACEOS file found')
  end if

  fileUnit = ut_getFreeFileUnit ()
  open (unit = fileUnit , file = tableName)
!
!
!   ...Read the temperature, density and energy group grids. Abort the calculation,
!      if any of the grids is not found. Check also, if the number of energy groups
!      read from the file corresponds to the one of the current run.
!
!
  do n = 1, 38
    read (fileUnit,*) 
  enddo

  nstepsTemperature = 0
  read (fileUnit,*) nstepsTemperature
  if (nstepsTemperature <= 0) then
      call Driver_abortFlash ('[op_readPropaceosTables] ERROR: no PROPACEOS temperature grid found')
  end if
  allocate(temperatures(nstepsTemperature))
  read (fileUnit,*) (temperatures(t),t=1,nstepsTemperature)

  nstepsDensity = 0
  read (fileUnit,*) nstepsDensity
  if (nstepsDensity <= 0) then
      call Driver_abortFlash ('[op_readPropaceosTables] ERROR: no PROPACEOS density grid found')
  end if
  allocate(densities(nstepsDensity))
  read (fileUnit,*) (densities(d),d=1,nstepsDensity)

! ... read and skip grid for opacities. It is the same as EOS grid
  do n = 1, 5
    read (fileUnit,*) 
  enddo
  read (fileUnit,*) nstepsTemperature
  read (fileUnit,*) (dummyData,t=1,nstepsTemperature)
  read (fileUnit,*) nstepsDensity
  read (fileUnit,*) (dummyData,d=1,nstepsDensity)

!
!
!   ...Read the energy group boundaries from the PROPACEOS file and compare with those stored.
!
!

  nEnergyGroups = 0
  read (fileUnit,*) nEnergyGroups
  if (nEnergyGroups <= 0) then
      call Logfile_stampMessage('[op_readPropaceosTables] WARNING: no PROPACEOS energy group grid found')
  end if
  read (fileUnit,*)
  read (fileUnit,*) (op_tabulatedEnergyBoundaries(g),g=1,nEnergyGroups+1)

  if (nEnergyGroups /= op_nEnergyGroups) then
     print *, nEnergyGroups, op_nEnergyGroups
      call Driver_abortFlash ('[op_readPropaceosTables] ERROR: bad size of PROPACEOS energy group grid')
  end if
  
  do g = 1,nEnergyGroups+1

     energyDifference = abs (op_tabulatedEnergyBoundaries (g) - op_energyGroupBoundaries (g)) / &
          op_energyGroupBoundaries(g)

     if (energyDifference > op_tableEnergyTolerance) then
         call Driver_abortFlash ('[op_readPropaceosTables] ERROR: PROPACEOS / Opacity energy group mismatch')
     end if
  end do
!
!
!   ...Perform the temperature eV -> K conversion.
!
!
  temperatures(:) = temperatures(:) * op_eV2Kelvin
!
!
!   ...Skip not needed data from the IONMIX4 file.
!
!
  notneededData =   nstepsDensity * nstepsTemperature

  read (fileUnit,*)
  read (fileUnit,*) (dummyData, n = 1,notneededData)
  read (fileUnit,*)
  read (fileUnit,*) (dummyData, n = 1,notneededData)
  read (fileUnit,*)
  read (fileUnit,*) (dummyData, n = 1,notneededData)
  read (fileUnit,*)
  read (fileUnit,*) (dummyData, n = 1,notneededData)
  read (fileUnit,*)
  read (fileUnit,*) (dummyData, n = 1,notneededData)
  read (fileUnit,*)
  read (fileUnit,*) (dummyData, n = 1,notneededData)
  read (fileUnit,*)
  read (fileUnit,*) (dummyData, n = 1,notneededData)
  read (fileUnit,*)
  read (fileUnit,*) (dummyData, n = 1,notneededData)
  read (fileUnit,*)
  read (fileUnit,*) (dummyData, n = 1,notneededData)
!
!
!   ...Establish the temperature and density grid for the Rosseland case (if needed)
!      and read in the opacities.
!
!
  do d=1, nstepsDensity
    do t=1, nstepsTemperature
      if (needROTable) then
        read (fileUnit,*)
        read (fileUnit,*) (op_RosselandTables (t,d,g,indexRO), g = 1,nEnergyGroups)
      else 
        read (fileUnit,*)
        read (fileUnit,*) (dummyData, g = 1,nEnergyGroups)
      end if
      if (needPETable) then
        read (fileUnit,*)
        read (fileUnit,*) (op_PlanckEmissionTables (t,d,g,indexPE), g = 1,nEnergyGroups)
      else 
        read (fileUnit,*)
        read (fileUnit,*) (dummyData, g = 1,nEnergyGroups)
      end if
      if (needPATable) then
        read (fileUnit,*)
        read (fileUnit,*) (op_PlanckAbsorptionTables (t,d,g,indexPA), g = 1,nEnergyGroups)
      else 
        read (fileUnit,*)
        read (fileUnit,*) (dummyData, g = 1,nEnergyGroups)
      end if
    enddo
  enddo
!
!
!   ...Establish the temperature and density grid for the Planck Absorption case (if needed)
!      and set  the opacities.
!
!
  if (needROTable) then

      op_nstepsDensityRO     (indexRO) = nstepsDensity
      op_nstepsTemperatureRO (indexRO) = nstepsTemperature

      op_tableDensityRO (1:nstepsDensity,indexRO) = densities(1:nstepsDensity)
      op_tableTemperatureRO (1:nstepsTemperature,indexRO) = temperatures(1:nstepsTemperature)

      if (op_useLogTables) then
                   op_tableDensityRO      (1:nstepsDensity,    indexRO) &
          = log10 (op_tableDensityRO      (1:nstepsDensity,    indexRO) )
                   op_tableTemperatureRO  (1:nstepsTemperature,indexRO) &
          = log10 (op_tableTemperatureRO  (1:nstepsTemperature,indexRO) )

          newZero = op_smallestPositiveNumber

                   op_RosselandTables     (1:nstepsTemperature,1:nstepsDensity,1:nEnergyGroups,indexRO) &
          =   max (op_RosselandTables     (1:nstepsTemperature,1:nstepsDensity,1:nEnergyGroups,indexRO),newZero)
                   op_RosselandTables     (1:nstepsTemperature,1:nstepsDensity,1:nEnergyGroups,indexRO) &
          = log10 (op_RosselandTables     (1:nstepsTemperature,1:nstepsDensity,1:nEnergyGroups,indexRO) )
      end if

  end if
!
!
!   ...Establish the temperature and density grid for the Planck Absorption case (if needed)
!      and read in the opacities.
!
!
  if (needPATable) then

      op_nstepsDensityPA     (indexPA) = nstepsDensity
      op_nstepsTemperaturePA (indexPA) = nstepsTemperature

      op_tableDensityPA (1:nstepsDensity,indexPA) = densities(1:nstepsDensity)
      op_tableTemperaturePA (1:nstepsTemperature,indexPA) = temperatures(1:nstepsTemperature)

      if (op_useLogTables) then
                   op_tableDensityPA         (1:nstepsDensity,    indexPA) &
          = log10 (op_tableDensityPA         (1:nstepsDensity,    indexPA) )
                   op_tableTemperaturePA     (1:nstepsTemperature,indexPA) &
          = log10 (op_tableTemperaturePA     (1:nstepsTemperature,indexPA) )

          newZero = op_smallestPositiveNumber

                   op_PlanckAbsorptionTables (1:nstepsTemperature,1:nstepsDensity,1:nEnergyGroups,indexPA) &
          =   max (op_PlanckAbsorptionTables (1:nstepsTemperature,1:nstepsDensity,1:nEnergyGroups,indexPA),newZero)
                   op_PlanckAbsorptionTables (1:nstepsTemperature,1:nstepsDensity,1:nEnergyGroups,indexPA) &
          = log10 (op_PlanckAbsorptionTables (1:nstepsTemperature,1:nstepsDensity,1:nEnergyGroups,indexPA) )
      end if

  end if
!
!
!   ...Establish the temperature and density grid for the Planck Emission case (if needed)
!      and read in the opacities.
!
!
  if (needPETable) then

      op_nstepsDensityPE     (indexPE) = nstepsDensity
      op_nstepsTemperaturePE (indexPE) = nstepsTemperature

      op_tableDensityPE (1:nstepsDensity,indexPE) = densities(1:nstepsDensity)
      op_tableTemperaturePE (1:nstepsTemperature,indexPE) = temperatures(1:nstepsTemperature)

      if (op_useLogTables) then
                   op_tableDensityPE       (1:nstepsDensity,    indexPE) &
          = log10 (op_tableDensityPE       (1:nstepsDensity,    indexPE) )
                   op_tableTemperaturePE   (1:nstepsTemperature,indexPE) &
          = log10 (op_tableTemperaturePE   (1:nstepsTemperature,indexPE) )

          newZero = op_smallestPositiveNumber

                   op_PlanckEmissionTables (1:nstepsTemperature,1:nstepsDensity,1:nEnergyGroups,indexPE) &
          =   max (op_PlanckEmissionTables (1:nstepsTemperature,1:nstepsDensity,1:nEnergyGroups,indexPE),newZero)
                   op_PlanckEmissionTables (1:nstepsTemperature,1:nstepsDensity,1:nEnergyGroups,indexPE) &
          = log10 (op_PlanckEmissionTables (1:nstepsTemperature,1:nstepsDensity,1:nEnergyGroups,indexPE) )
      end if

  end if

!
!
!   ...Close the PROPACEOS file and deallocate temperature/density arrays
!
!
  deallocate(densities)
  deallocate(temperatures)
  close (fileUnit)
!
!
!   ...Ready! 
!
!
  return
end subroutine op_readPropaceosTables
