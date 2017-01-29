!!****if* source/Simulation/SimulationMain/unitTest/SolidBody/MultipleWithTriangles/MovingBodies/NoComm/gr_ptAdvance
!!
!! NAME
!!
!!  gr_ptAdvance
!!
!! SYNOPSIS
!!
!!  call gr_ptAdvance(real(in)   :: dtOld,
!!                         real(in)   :: dtNew)
!!
!! DESCRIPTION
!!
!!  Time advancement routine for the body.
!!  
!! ARGUMENTS
!!
!!   dtOld -- not used in this first-order scheme
!!   dtNew -- current time increment

#include "constants.h"
#include "Flash.h"

subroutine gr_ptAdvance (dtOld,dtNew)

  use gr_sbData, ONLY : gr_sbBodyInfo, NumVertices, gr_sbNumBodies
  implicit none

  real, INTENT(in)  :: dtOld, dtNew
  integer       :: i, b
  real          :: jumpx,jumpy,jumpz

  ! Update the particle positions.
  do b = 1, gr_sbNumBodies
     jumpx = dtNew 
     jumpy = dtNew 
     jumpz = dtNew 
     do i = 1, NumVertices
        gr_sbBodyInfo(b) % xb(i) = gr_sbBodyInfo(b) % xb(i) + jumpx
        gr_sbBodyInfo(b) % yb(i) = gr_sbBodyInfo(b) % yb(i) + jumpy
        gr_sbBodyInfo(b) % zb(i) = gr_sbBodyInfo(b) % zb(i) + jumpz
     enddo
  enddo
  return
end subroutine gr_ptAdvance


