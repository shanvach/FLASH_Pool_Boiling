!!****if* source/Simulation/SimulationMain/INavierStokes/2D/IB_Cyl_parallelIBVP/old/ImBound
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
  use Grid_interface, only : Grid_getBoundboxCentroids,Grid_sbSelectMaster
  use gr_sbInterface, only : gr_sbCreateParticles,gr_sbGetProcBlock, &
                             gr_sbSendPosn, gr_sbDistributedForces
  use ib_interface, only : ib_forcing,ib_CalcForce

  use ImBound_data, only : ib_dt

  use gr_sbData, only : gr_sbNumBodies,gr_sbBodyInfo

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


  ib_dt = dt

  select case (forcflag)
  case(FORCE_FLOW) 

  ! First, select Master:
  call Grid_getBoundboxCentroids()
  call Grid_sbSelectMaster()
  call gr_sbCreateParticles()
  call gr_sbGetProcBlock()
  call gr_sbSendPosn()

  case(COMPUTE_FORCES)

  call gr_sbDistributedForces()
  
  case default

  call Driver_abortFlash("ImBound : forcflag doen not correspond to any available option.") 
  
  end select

  return
end subroutine ImBound

