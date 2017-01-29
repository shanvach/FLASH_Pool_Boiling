!!****if* source/physics/sourceTerms/Cool/CoolMain/Radloss/Cool_init
!!
!! NAME
!!
!!  Cool_init
!!
!! SYNOPSIS
!!
!!  Cool_init()
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   
!!
!!
!!***

subroutine Cool_init()
  use Cool_data, only : cl_meshMe, cl_Xin, cl_Abar,&
       cl_tradmin, cl_tradmax, cl_dradmin, cl_dradmax,useCool, &
       cl_h1SpecIndex, cl_elecSpecIndex, &
       cl_speciesNameH1, cl_speciesNameElec
  use Simulation_interface, ONLY : Simulation_mapStrToInt
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Multispecies_interface, ONLY : Multispecies_getProperty
  use Grid_interface, ONLY: Grid_getNumProcs

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"

  implicit none
  
  
  integer :: i

  ! Everybody should know these
  call Driver_getMype(MESH_COMM,cl_meshMe)

  call RuntimeParameters_get("tradmin", cl_tradmin)
  call RuntimeParameters_get("tradmax", cl_tradmax)
  call RuntimeParameters_get("dradmin", cl_dradmin)
  call RuntimeParameters_get("dradmax", cl_dradmax)
  call RuntimeParameters_get("useCool",useCool)
  if (.not. useCool) then
     write(6,*)'WARNING:  You have included the Cool unit but have set '
     write(6,*)'   the runtime parameter useCool to FALSE'
     write(6,*)'   No cooling will occur but Cool_init will continue.'
  end if

  call RuntimeParameters_get("cl_speciesNameH1",   cl_speciesNameH1)
  call RuntimeParameters_get("cl_speciesNameElec", cl_speciesNameElec)
  call Simulation_mapStrToInt(cl_speciesNameH1,  cl_h1SpecIndex,MAPBLOCK_UNK)
  call Simulation_mapStrToInt(cl_speciesNameElec,cl_elecSpecIndex,MAPBLOCK_UNK)

  print*,'cl_h1SpecIndex, cl_elecSpecIndex:',cl_h1SpecIndex, cl_elecSpecIndex

  do i = 1, NSPECIES
     call Multispecies_getProperty(SPECIES_BEGIN+i-1,A, cl_Abar(i))
  end do
!        make sure composition has protons and electrons


!===========================================================================

  return
end subroutine Cool_init
