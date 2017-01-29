!!****if* source/Simulation/SimulationMain/INavierStokes/2D/IB_Cyl_parallelIB/ImBound
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
                           ib_ubd,ib_vbd,ib_ubdd,ib_vbdd,omg,freq_t,Ro,ao,xo,yo,ib_dt

  use gr_sbData, only : gr_sbNumBodies,gr_sbBodyInfo

  use Grid_interface, only : Grid_sbSelectMaster
  use gr_sbInterface, only : gr_sbCreateParticles,gr_sbGetProcBlock, &
                             gr_sbSendPosn

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

  ib_dt = dt

  select case (forcflag)
  case(FORCE_FLOW) 

  ! Do surf kinematics gr_ptAdvance():
  tita  = omg*dr_simTime
  do ibd=1,gr_sbNumBodies
  do i=1,gr_sbBodyInfo(ibd)%NumVertices


    ! Positions of markers, normals and arclength var
    gr_sbBodyInfo(ibd)%xb(i) =  gr_sbBodyInfo(ibd)%xo +  gr_sbBodyInfo(ibd)%xbus(i)*cos(tita) &
                               -gr_sbBodyInfo(ibd)%ybus(i)*sin(tita)

    
    gr_sbBodyInfo(ibd)%yb(i) =  gr_sbBodyInfo(ibd)%yo +  gr_sbBodyInfo(ibd)%xbus(i)*sin(tita) &
                              + gr_sbBodyInfo(ibd)%ybus(i)*cos(tita) +ao*sin(2.*PI*freq_t*dr_simTime)

    
    gr_sbBodyInfo(ibd)%nxL(i) =  gr_sbBodyInfo(ibd)%xb(i)/Ro
    gr_sbBodyInfo(ibd)%nyL(i) =  gr_sbBodyInfo(ibd)%yb(i)/Ro

    ! Set velocities and accelerations
    gr_sbBodyInfo(ibd)%ubd(i) = -omg*gr_sbBodyInfo(ibd)%xbus(i)*sin(tita) - &
                                 omg*gr_sbBodyInfo(ibd)%ybus(i)*cos(tita) ! 0.
    gr_sbBodyInfo(ibd)%vbd(i) =  omg*gr_sbBodyInfo(ibd)%xbus(i)*cos(tita) - &
                                 omg*gr_sbBodyInfo(ibd)%ybus(i)*sin(tita) + &
                                 (2.*PI*freq_t)*ao*cos(2.*PI*freq_t*dr_simTime) ! 0.
    gr_sbBodyInfo(ibd)%ubdd(i)= -omg**2.*gr_sbBodyInfo(ibd)%xbus(i)*cos(tita) + &
                                 omg**2.*gr_sbBodyInfo(ibd)%ybus(i)*sin(tita)
    gr_sbBodyInfo(ibd)%vbdd(i)= -omg**2.*gr_sbBodyInfo(ibd)%xbus(i)*sin(tita) - &
                                 omg**2.*gr_sbBodyInfo(ibd)%ybus(i)*cos(tita) - &
                                 (2.*PI*freq_t)**2.*ao*sin(2.*PI*freq_t*dr_simTime)



  enddo 
  enddo


!  call Grid_getBoundboxCentroids()
  call Grid_sbSelectMaster()
  call gr_sbCreateParticles()
  call gr_sbGetProcBlock()
  call gr_sbSendPosn()

  !call ib_forcing(blockCount, blockList, dt)

  case(COMPUTE_FORCES)

  !call ib_CalcForce(blockCount, blockList, dt)

  case default

  call Driver_abortFlash("ImBound : forcflag doen not correspond to any available option.") 
  
  end select

  return
end subroutine ImBound

