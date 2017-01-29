!! NAME
!!
!!
!!
!! SYNOPSIS
!!  
!! VARIABLES
!!
!!
!! DESCRIPTION
!! Set the parameters for the collision model 
!! shizhao Wang
!! Nov 05, 2014
!!***

#include "constants.h"
#include "Flash.h"
#include "SolidMechanics.h"

subroutine sm_set_cm_para()

  use SolidMechanics_data, only : sm_meshMe, sm_BodyInfo, sm_cm_rhou, sm_cm_eps1, sm_cm_eps2
  use Grid_interface, ONLY : Grid_getMinCellSize
  use Driver_interface, only : Driver_abortFlash

  implicit none

  real :: Dmin

#ifdef FLASH_GRID_PARAMESH 
  call Driver_abortFlash('The routine sm_set_cm_para has not been tested for the PARAMESH')  
#else /* Uniform Grid */
  call Grid_getMinCellSize(Dmin)
#endif

  sm_cm_rhou = 2.*Dmin  ! because the exteral probe point is at nh=2*Dmin
  sm_cm_eps1 = Dmin*Dmin/100.
  sm_cm_eps2 = Dmin*Dmin/100.

  if(sm_meshMe == MASTER_PE) then
    write(*,*) 'Parameters for the collision model'
    write(*,*) 'minimum grid size:', Dmin
    write(*,*) 'sm_cm_rhou:', sm_cm_rhou
    write(*,*) 'sm_cm_eps1:', sm_cm_eps1
    write(*,*) 'sm_cm_eps2:', sm_cm_eps2
  endif
  
  return

end subroutine sm_set_cm_para
