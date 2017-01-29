!!****if* source/Simulation/SimulationMain/unitTest/SolidBody/SingleWithTriangles/MovingBody/Grid_getBoundboxCentroids
!!
!! NAME
!!  Grid_getBoundboxCentroids
!!
!! SYNOPSIS
!!
!!  Grid_getBoundboxCentroids()
!!
!! DESCRIPTION
!!
!!  Called from Driver_evolveFlash for moving body. Calculate the boundbox of each solid body.  
!!  Also calculate the centriod of each triangle within each solid body
!!
!! ARGUMENTS

#include "constants.h"
#include "Flash.h"

Subroutine Grid_getBoundboxCentroids()
  use Grid_data, ONLY : gr_meshMe
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, &
                        gr_sbDebug, &
                        aelem, NumTriangles

  implicit none
  integer :: i

  gr_sbBodyInfo(1) % boundBox(LOW,IAXIS) = minval(gr_sbBodyInfo(1) % xb)
  gr_sbBodyInfo(1) % boundBox(HIGH,IAXIS) = maxval(gr_sbBodyInfo(1) % xb)
  gr_sbBodyInfo(1) % boundBox(LOW,JAXIS) = minval(gr_sbBodyInfo(1) % yb)
  gr_sbBodyInfo(1) % boundBox(HIGH,JAXIS) = maxval(gr_sbBodyInfo(1) % yb)
  gr_sbBodyInfo(1) % boundBox(LOW,KAXIS) = minval(gr_sbBodyInfo(1) % zb)
  gr_sbBodyInfo(1) % boundBox(HIGH,KAXIS) = maxval(gr_sbBodyInfo(1) % zb)

  !Calculating centroids for each triangle
  do i = 1, NumTriangles
     gr_sbBodyInfo(1) % triangleCentroids(IAXIS,i) = &
          1./3. * (gr_sbBodyInfo(1) % xb(aelem(1,i)) + gr_sbBodyInfo(1) % xb(aelem(2,i)) + gr_sbBodyInfo(1) % xb(aelem(3,i)))

     gr_sbBodyInfo(1) % triangleCentroids(JAXIS,i) = &
          1./3. * (gr_sbBodyInfo(1) % yb(aelem(1,i)) + gr_sbBodyInfo(1) % yb(aelem(2,i)) + gr_sbBodyInfo(1) % yb(aelem(3,i)))

     gr_sbBodyInfo(1) % triangleCentroids(KAXIS,i) = &
          1./3. * (gr_sbBodyInfo(1) % zb(aelem(1,i)) + gr_sbBodyInfo(1) % zb(aelem(2,i)) + gr_sbBodyInfo(1) % zb(aelem(3,i)))
  end do
  return
End Subroutine Grid_getBoundboxCentroids
