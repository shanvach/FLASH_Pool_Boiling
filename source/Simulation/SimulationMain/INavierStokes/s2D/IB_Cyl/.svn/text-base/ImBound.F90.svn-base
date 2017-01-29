!!****if* source/Simulation/SimulationMain/INavierStokes/2D/IB_Cyl/ImBound
!!
!!
!! NAME
!!
!!  ImBound
!!
!!
!! SYNOPSIS
!!
!!  ImBound(blockCount,blockList,dt)
!!
!!
!! DESCRIPTION
!!
!!
!!
!!***

subroutine ImBound(blockCount, blockList, dt, forcflag)

  use Driver_data, only : dr_simTime
  use Driver_interface, only : Driver_abortFlash

  use ib_interface, only : ib_forcing,ib_CalcForce

  use ImBound_data, only : ib_ABODY,ib_nbd,ib_xb,ib_yb,ib_xbus,ib_ybus,ib_nxL,ib_nyL, &
                           ib_ubd,ib_vbd,ib_ubdd,ib_vbdd,omg,freq_t,Ro,ao,xo,yo
  implicit none
#include "Flash.h"
#include "ImBound.h"
#include "constants.h"

  !! ---- Argument List ----------------------------------
  integer, INTENT(IN) :: blockCount
  integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList
  real, INTENT(IN) :: dt
  integer, INTENT(IN) :: forcflag
  !! -----------------------------------------------------

  real :: tita
  integer ibd,i

  select case (forcflag)
  case(FORCE_FLOW) 

  ! Do surf kinematics:
  tita  = omg*dr_simTime
  do ibd=1,ib_nbd
  do i=ib_ABODY(ibd)%lb+1,ib_ABODY(ibd)%lb+ib_ABODY(ibd)%mb


    ! Positions of markers, normals and arclength var
    ib_xb(i) = xo(ibd) + ib_xbus(i)*cos(tita) - ib_ybus(i)*sin(tita)
    ib_yb(i) = yo(ibd) + ib_xbus(i)*sin(tita) + ib_ybus(i)*cos(tita) +ao*sin(2.*PI*freq_t*dr_simTime)

    ib_nxL(i) = ib_xb(i)/Ro
    ib_nyL(i) = ib_yb(i)/Ro

    ! Set velocities and accelerations
    ib_ubd(i) = -omg*ib_xbus(i)*sin(tita) - omg*ib_ybus(i)*cos(tita) ! 0.
    ib_vbd(i) =  omg*ib_xbus(i)*cos(tita) - omg*ib_ybus(i)*sin(tita) +(2.*PI*freq_t)*ao*cos(2.*PI*freq_t*dr_simTime) ! 0.
    ib_ubdd(i)= -omg**2.*ib_xbus(i)*cos(tita) + omg**2.*ib_ybus(i)*sin(tita)
    ib_vbdd(i)= -omg**2.*ib_xbus(i)*sin(tita) - omg**2.*ib_ybus(i)*cos(tita) - (2.*PI*freq_t)**2.*ao*sin(2.*PI*freq_t*dr_simTime)

  enddo 
  enddo


  call ib_forcing(blockCount, blockList, dt)

  case(COMPUTE_FORCES)

  call ib_CalcForce(blockCount, blockList, dt)

  case default

  call Driver_abortFlash("ImBound : forcflag doen not correspond to any available option.") 
  
  end select

  return
end subroutine ImBound

