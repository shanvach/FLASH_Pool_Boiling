!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Tabulated/op_browsePropaceosTables
!!
!! NAME
!!
!!  op_browsePropaceosTables
!!
!! SYNOPSIS
!!
!!  call op_browsePropaceosTables (character (in)  :: tableName (len=80),
!!                                 logical   (in)  :: needPATable,
!!                                 logical   (in)  :: needPETable,
!!                                 logical   (in)  :: needROTable,
!!                                 integer   (out) :: nstepsDensityPA,
!!                                 integer   (out) :: nstepsDensityPE,
!!                                 integer   (out) :: nstepsDensityRO,
!!                                 integer   (out) :: nstepsTemperaturePA,
!!                                 integer   (out) :: nstepsTemperaturePE,
!!                                 integer   (out) :: nstepsTemperatureRO)
!!
!! DESCRIPTION
!!
!!  This routine browses through the tabulated opacities from an PROPACEOS datafile output in
!!  order to extract the number of steps for both the density and the temperature grid with
!!  which the PROPACEOS tables were generated.
!!
!! ARGUMENTS
!!
!!  tableName           : the name of the PROPACEOS file
!!  needPATable         : if yes, Planck Absorption Opacities are needed from the PROPACEOS table
!!  needPETable         : if yes, Planck   Emission Opacities are needed from the PROPACEOS table
!!  needROTable         : if yes,         Rosseland Opacities are needed from the PROPACEOS table
!!  nstepsDensityPA     : the size of the Planck Absorption density grid returned
!!  nstepsDensityPE     : the size of the Planck   Emission density grid returned
!!  nstepsDensityRO     : the size of the         Rosseland density grid returned
!!  nstepsTemperaturePA : the size of the Planck Absorption temperature grid returned
!!  nstepsTemperaturePE : the size of the Planck   Emission temperature grid returned
!!  nstepsTemperatureRO : the size of the         Rosseland temperature grid returned
!!
!!***

#include "constants.h"

subroutine op_browsePropaceosTables (tableName,                        &
                                     needPATable,                      &
                                     needPETable,                      &
                                     needROTable,                      &
                                                  nstepsDensityPA,     &
                                                  nstepsDensityPE,     &
                                                  nstepsDensityRO,     &
                                                  nstepsTemperaturePA, &
                                                  nstepsTemperaturePE, &
                                                  nstepsTemperatureRO  )

  use Driver_interface,  ONLY : Driver_abortFlash
  use Opacity_data,  ONLY : op_globalMe

  implicit none

  character (len=80), intent (in)  :: tableName
  logical,            intent (in)  :: needPATable
  logical,            intent (in)  :: needPETable
  logical,            intent (in)  :: needROTable
  integer,            intent (out) :: nstepsDensityPA
  integer,            intent (out) :: nstepsDensityPE
  integer,            intent (out) :: nstepsDensityRO
  integer,            intent (out) :: nstepsTemperaturePA
  integer,            intent (out) :: nstepsTemperaturePE
  integer,            intent (out) :: nstepsTemperatureRO

  logical :: fileExists

  integer :: fileUnit,  i
  integer :: nstepsDensity
  integer :: nstepsTemperature
  integer :: ut_getFreeFileUnit
  
  real    :: dummy
!
!
!   ...Check and open the PROPACEOS opacity file.
!
!
  inquire (file = tableName , exist = fileExists)

  if (.not.fileExists) then
     if (op_globalMe==MASTER_PE) &
          print*,'[op_browsePropaceosTables] ERROR: PROPACEOS file not found: ',tableName 
      call Driver_abortFlash ('[op_browsePropaceosTables] ERROR: no PROPACEOS file found')
  end if

  fileUnit = ut_getFreeFileUnit ()
  open (unit = fileUnit , file = tableName)
!
!
!   ...Read the temperature and density grids. Abort the calculation,
!      if any of the grids is not found.
!
!
  do i = 1, 38
    read (fileUnit,*) 
  enddo
  
  nstepsTemperature = 0
  read (fileUnit,*) nstepsTemperature
  if (nstepsTemperature <= 0) then
      call Driver_abortFlash ('[eos_tabBrowsePropaceosTables] ERROR: no PROPACEOS temperature grid found')
  end if
  read (fileUnit,*) (dummy,i=1,nstepsTemperature)

  nstepsDensity = 0
  read (fileUnit,*) nstepsDensity
  if (nstepsDensity <= 0) then
      call Driver_abortFlash ('[eos_tabBrowsePropaceosTables] ERROR: no PROPACEOS density grid found')
  end if
  read (fileUnit,*) (dummy,i=1,nstepsDensity)

!
!
!   ...For the PROPACEOS tables the grids are the same for all three Planck Absorption,
!      Planck Emission and Rosseland Transporation tables.
!
!
  nstepsDensityPA = nstepsDensity
  nstepsDensityPE = nstepsDensity
  nstepsDensityRO = nstepsDensity

  nstepsTemperaturePA = nstepsTemperature
  nstepsTemperaturePE = nstepsTemperature
  nstepsTemperatureRO = nstepsTemperature
!
!
!   ...Close the PROPACEOS file.
!
!
  close (fileUnit)
!
!
!   ...Ready! 
!
!
  return
end subroutine op_browsePropaceosTables
