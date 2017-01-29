!!****if* source/physics/sourceTerms/Cool/CoolMain/Isothermal/Cool_init
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
!! AUTOGENROBODOC
!!
!!
!!***

subroutine Cool_init()
  use Cool_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
#include "constants.h"
#include "Eos.h"
  implicit none
  
  integer :: vecLen
 
  
  call RuntimeParameters_get ("isotherm", cl_isotherm)
  call RuntimeParameters_get ("smallt", cl_smallt)
  call RuntimeParameters_get ("smallp", cl_smallp)
  call RuntimeParameters_get ("smalle", cl_smalle)
  call RuntimeParameters_get ("useCool",useCool)
  if (.not. useCool) then
     write(6,*)'WARNING:  You have included the Cool unit but have set '
     write(6,*)'   the runtime parameter useCool to FALSE'
     write(6,*)'   No cooling will occur but Cool_init will continue.'
  end if
  call Driver_getMype(MESH_COMM,cl_meshMe)

  return
end subroutine Cool_init
