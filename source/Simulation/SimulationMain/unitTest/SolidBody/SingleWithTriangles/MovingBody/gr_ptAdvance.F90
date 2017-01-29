!!****if* source/Simulation/SimulationMain/unitTest/SolidBody/SingleWithTriangles/MovingBody/gr_ptAdvance
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

  use gr_sbData, ONLY : gr_sbBodyInfo, NumVertices
  implicit none

  real, INTENT(in)  :: dtOld, dtNew
  integer       :: i
  real          :: jumpx,jumpy,jumpz

  ! Update the particle positions.
  jumpx = dtNew * gr_sbBodyInfo(1) % boundBox(HIGH,IAXIS) / 2
  jumpy = dtNew * gr_sbBodyInfo(1) % boundBox(HIGH,JAXIS) / 2
  jumpz = dtNew * gr_sbBodyInfo(1) % boundBox(HIGH,KAXIS) / 2
  do i = 1, NumVertices
     gr_sbBodyInfo(1) % xb(i) = gr_sbBodyInfo(1) % xb(i) + jumpx
     gr_sbBodyInfo(1) % yb(i) = gr_sbBodyInfo(1) % yb(i) + jumpy
     gr_sbBodyInfo(1) % zb(i) = gr_sbBodyInfo(1) % zb(i) + jumpz
  enddo
  return
end subroutine gr_ptAdvance


