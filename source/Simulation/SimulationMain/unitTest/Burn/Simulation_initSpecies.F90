!!****if* source/Simulation/SimulationMain/unitTest/Burn/Simulation_initSpecies
!!
!! NAME
!!
!!  Simulation_initSpecies
!!
!!
!! SYNOPSIS
!!  Simulation_initSpecies()
!!
!! DESCRIPTION
!!
!!  This routine will initialize the species and species values needed for
!!  the burning unit test.  
!!
!! NOTES
!!  The two species names are defined as 
!!  PARAMETER sim_compA and sim_compB in the Config file and then given in the
!!  flash.par.  This routine translates the string name into a proper species,
!!  similar to an inverse of SimulationComposition.  New indices are stored as
!!  sim_compIndexA and B
!!
!!  Unbeknownst to LBR, Simulation_initSpecies is executed BEFORE Simulation_init.
!!  So you'll have to do RuntimeParameter initialization here.....
!!
!!***

subroutine Simulation_initSpecies()

  use Simulation_data, ONLY: sim_compA, sim_compB, sim_compIndexA, sim_compIndexB
  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Simulation_interface, ONLY : Simulation_mapStrToInt
  use Multispecies_interface, ONLY : Multispecies_setProperty

  implicit none

#include "constants.h"
#include "Multispecies.h"
#include "Flash.h"

  integer, parameter :: SPEC_UNIT=2,SPEC_NUM=238
  character(len=4) :: isotopeName
  integer :: isotopeInteger, i
  real :: abar, zbar, bindEnergy, dummy
  logical :: foundA, foundB
  character*90 :: internalFile

  !! Note that Simulation_initSpecies is called BEFORE Simulation_init
  call RuntimeParameters_get("compA", sim_compA)
  call RuntimeParameters_get("compB", sim_compB)


  !! Convert the user's request into integers
  call Simulation_mapStrToInt(sim_compA,sim_compIndexA,MAPBLOCK_UNK)
  if (sim_compIndexA == NONEXISTENT) then
     write(internalFile,*)'Could not find species',sim_compA,' in Config variables'
     print *, internalFile
     call Driver_abortFlash(internalFile)
  endif
  foundA = .false.          !! still need to map it to something from SpeciesList.txt file

  call Simulation_mapStrToInt(sim_compB,sim_compIndexB,MAPBLOCK_UNK)
  if (sim_compIndexB == NONEXISTENT) then
     write(internalFile,*)'Could not find species',sim_compB,' in Config variables'
     print *, internalFile
     call Driver_abortFlash(internalFile)
  endif
  foundB = .false.

  open(unit=SPEC_UNIT,file="SpeciesList.txt")

  do i=1, SPEC_NUM
     read(2,*)isotopeName,dummy,abar,zbar,bindEnergy
     call Simulation_mapStrToInt(isotopeName,isotopeInteger,MAPBLOCK_UNK)
     if(isotopeInteger /= NONEXISTENT) then    !! must be either compa or compB
        if (isotopeInteger == sim_compIndexA)  foundA = .true.
        if (isotopeInteger == sim_compIndexB)  foundB = .true.
        call Multispecies_setProperty(isotopeInteger, A, abar)
        call Multispecies_setProperty(isotopeInteger, Z, zbar)
        call Multispecies_setProperty(isotopeInteger, EB, bindEnergy)
     end if
  end do

  if ((.NOT. foundA) .AND. (.NOT.foundB)) then
     write(internalFile,*)'Could not find isotope information for ',sim_compA,' or ',sim_compB
     call Driver_abortFlash(internalFile)
  endif

  return

end subroutine Simulation_initSpecies
