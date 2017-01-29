!!****if* source/physics/sourceTerms/Cool/CoolMain/SutherlandDopita/Cool_computeDt
!!
!! NAME
!!  
!!  Cool_computeDt
!!
!!
!! SYNOPSIS
!! 
!!  Cool_computeDt ( integer(IN) : blockID, 
!!                   
!!                  real,pointer :  solnData(:,:,:,:),   
!!                  real,(INOUT):   dt_check, 
!!                  integer(INOUT): dt_minloc(:) )
!!  
!! DESCRIPTION
!!
!!  Computes the timestep limiter for cooling source term solver.
!!  Version for Sutherland & Dopita tabulated cooling function
!!
!!
!! ARGUMENTS
!!
!!  blockID        local block ID
!!  
!!  solnData        the physical, solution data from grid
!!  dt_check        variable to hold timestep constraint
!!  dt_minloc(5)    array to hold limiting zone info:  zone indices
!!
!!***

subroutine Cool_computeDt (blockID, &
                              blkLimits,blkLimitsGC,        &
                              solnData,   &
                              dt_check, dt_minloc )
  use Driver_interface, ONLY : Driver_abortFlash
  use Cool_data, ONLY : useCool,cl_rho, cl_gamma, cl_meshMe
#include "Flash.h"
#include "constants.h"

  implicit none

  integer, intent(IN) :: blockID
  integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
  real,INTENT(INOUT)    :: dt_check
  integer,INTENT(INOUT)    :: dt_minloc(5)
  real, pointer, dimension(:,:,:,:) :: solnData
  

  real    :: dedt, eid, dt_temp, tcool
  integer :: i, j, k, temploc(5)
  
  !==============================================================================
  

  if(.not.useCool) return

! Initialize dt_temp to some obscenely large number
  dt_temp = HUGE(1.0)

! sweep over all the zones
  
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
           
           cl_rho = solnData(DENS_VAR,i,j,k)
           cl_gamma = solnData(GAME_VAR,i,j,k)
           eid = solnData(DENS_VAR,i,j,k) * solnData(EINT_VAR,i,j,k)
   
           call cool_deriv (eid, dedt)
           
           ! use 1/2 the cooling time as the limiter
           tcool = 0.5 * eid / max(abs(dedt), 1.E-99)

           if (tcool < dt_temp) then
              dt_temp = tcool
              temploc(1) = i
              temploc(2) = j
              temploc(3) = k
              temploc(4) = blockID
              temploc(5) = cl_meshMe
           endif
           
        enddo
     enddo
  enddo
  
  if (dt_temp < dt_check) then
     dt_check = dt_temp
     dt_minloc = temploc
  endif
  
  if(dt_check <= 0.0) call Driver_abortFlash("[Cool]: computed dt is not positive! Aborting!")
!==============================================================================

  return
end subroutine Cool_computeDt

