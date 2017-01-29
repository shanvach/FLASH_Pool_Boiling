! sm_pc_EstDT_3DFlexible.F90
!
! Search over each element of a 3DFlexible body and estimate the DT 
! by approximating the largest eigenvalue of each element.

#include "Flash.h"
#include "constants.h"
#include "SolidMechanics.h"

subroutine sm_EstDT_3DFlexible(ibd, dt)

  use SolidMechanics_data,  only: sm_structure, sm_BodyInfo
  use sm_element_interface, only: sm_3DFlexible_getElement_EvalMax
  use Driver_interface,     only: Driver_abortFlash
  implicit none

  ! IO Variables
  integer, intent(in)  :: ibd
  real,    intent(out) :: dt
  
  ! Hack
  dt = 1.e-3

end subroutine sm_EstDT_3DFlexible

