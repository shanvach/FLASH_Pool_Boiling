!!****if* source/Simulation/SimulationMain/unitTest/SolidBody/MultipleWithTriangles/MovingBodies/Grid_getBoundboxCentroids
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
  integer :: b, i

  do b = 1, gr_sbNumBodies
     gr_sbBodyInfo(b) % boundBox(LOW,IAXIS) = minval(gr_sbBodyInfo(b) % xb)
     gr_sbBodyInfo(b) % boundBox(HIGH,IAXIS) = maxval(gr_sbBodyInfo(b) % xb)
     gr_sbBodyInfo(b) % boundBox(LOW,JAXIS) = minval(gr_sbBodyInfo(b) % yb)
     gr_sbBodyInfo(b) % boundBox(HIGH,JAXIS) = maxval(gr_sbBodyInfo(b) % yb)
     gr_sbBodyInfo(b) % boundBox(LOW,KAXIS) = minval(gr_sbBodyInfo(b) % zb)
     gr_sbBodyInfo(b) % boundBox(HIGH,KAXIS) = maxval(gr_sbBodyInfo(b) % zb)
  
     !Calculating centroids for each triangle
!     do i = 1, NumTriangles
!        gr_sbBodyInfo(b) % triangleCentroids(IAXIS,i) = &
!             1./3. * (gr_sbBodyInfo(b) % xb(aelem(1,i)) + gr_sbBodyInfo(b) % xb(aelem(2,i)) + gr_sbBodyInfo(b) % xb(aelem(3,i)))
        
!        gr_sbBodyInfo(b) % triangleCentroids(JAXIS,i) = &
!             1./3. * (gr_sbBodyInfo(b) % yb(aelem(1,i)) + gr_sbBodyInfo(b) % yb(aelem(2,i)) + gr_sbBodyInfo(b) % yb(aelem(3,i)))
     
!        gr_sbBodyInfo(b) % triangleCentroids(KAXIS,i) = &
!             1./3. * (gr_sbBodyInfo(b) % zb(aelem(1,i)) + gr_sbBodyInfo(b) % zb(aelem(2,i)) + gr_sbBodyInfo(b) % zb(aelem(3,i)))
!     end do
  enddo
  return
End Subroutine Grid_getBoundboxCentroids
