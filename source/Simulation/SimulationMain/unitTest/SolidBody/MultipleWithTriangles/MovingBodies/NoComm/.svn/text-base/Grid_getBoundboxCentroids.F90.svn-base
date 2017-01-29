!!****if* source/Simulation/SimulationMain/unitTest/SolidBody/MultipleWithTriangles/MovingBodies/NoComm/Grid_getBoundboxCentroids
!!
!! NAME
!!  Grid_getBoundboxCentroids
!!
!! SYNOPSIS
!!
!!  Grid_getBoundboxCentroids
!!
!! DESCRIPTION
!!
!!  Called from gr_sbInit. Calculate the boundbox of each solid body.  
!!
!! ARGUMENTS

#include "constants.h"
#include "Flash.h"

Subroutine Grid_getBoundboxCentroids()
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies
 
  implicit none
  integer :: b

  do b = 1, gr_sbNumBodies
     gr_sbBodyInfo(b) % boundBox(LOW,IAXIS) = minval(gr_sbBodyInfo(b) % xb)
     gr_sbBodyInfo(b) % boundBox(HIGH,IAXIS) = maxval(gr_sbBodyInfo(b) % xb)
     gr_sbBodyInfo(b) % boundBox(LOW,JAXIS) = minval(gr_sbBodyInfo(b) % yb)
     gr_sbBodyInfo(b) % boundBox(HIGH,JAXIS) = maxval(gr_sbBodyInfo(b) % yb)
     gr_sbBodyInfo(b) % boundBox(LOW,KAXIS) = minval(gr_sbBodyInfo(b) % zb)
     gr_sbBodyInfo(b) % boundBox(HIGH,KAXIS) = maxval(gr_sbBodyInfo(b) % zb)

  enddo
  return
End Subroutine Grid_getBoundboxCentroids
