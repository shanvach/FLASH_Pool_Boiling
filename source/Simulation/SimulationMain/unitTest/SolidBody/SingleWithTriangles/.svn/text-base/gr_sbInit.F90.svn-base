!!****if* source/Simulation/SimulationMain/unitTest/SolidBody/SingleWithTriangles/gr_sbInit
!!
!! NAME
!!  gr_sbInit
!!
!! SYNOPSIS
!!
!!  gr_sbInit()
!!
!! DESCRIPTION
!!
!!  * Called from Grid_init
!!  * Read input file containing the coordinates of the vertices and the triangle elements.
!!  * Get the boundary box of the body
!!  * Create triangles representing the body
!!  * Calculate centroids of each triangle
!!
!! ARGUMENTS
!! 

#include "constants.h"
#include "Flash.h"

Subroutine gr_sbInit()
  use Grid_data, ONLY : gr_meshMe
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, &
       gr_sbPtNumX, gr_sbPtNumY, gr_sbPtNumZ, gr_sbDebug, &
       aelem, NumTriangles, NumVertices
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none
  integer :: i
  real, allocatable, dimension(:) :: xb_elem, yb_elem, zb_elem

  open (unit=9, file = 'RBC.dat', status = 'old', action = 'read')
  read (9, '(I4)') NumVertices

  read (9, '(I4)') NumTriangles

  gr_sbNumBodies = 1
  allocate(gr_sbBodyInfo(gr_sbNumBodies))

  allocate(gr_sbBodyInfo(1) % xb(NumVertices), gr_sbBodyInfo(1) % yb(NumVertices), gr_sbBodyInfo(1) % zb(NumVertices)) !Vertex points
  allocate(aelem(3,NumTriangles))  !Triangles
  allocate(xb_elem(NumTriangles), yb_elem(NumTriangles), zb_elem(NumTriangles)) !Vertex numbers of trianlge

  gr_sbBodyInfo(1) % boundBox(:,:) = 0.0

  call RuntimeParameters_get("sb_ptNumX", gr_sbPtNumX)

  if (NDIM >= 2) then
     call RuntimeParameters_get("sb_ptNumY", gr_sbPtNumY)
  else
     gr_sbPtNumY = 1
  end if

  if (NDIM == 3) then
     call RuntimeParameters_get("sb_ptNumZ", gr_sbPtNumZ)
  else
     gr_sbPtNumZ = 1
  end if

  do i = 1, NumVertices
     read (9, *) gr_sbBodyInfo(1) % xb(i), gr_sbBodyInfo(1) % yb(i), gr_sbBodyInfo(1) % zb(i) !read position coordinates of vertices
  enddo
  gr_sbBodyInfo(1) % boundBox(LOW,IAXIS) = minval(gr_sbBodyInfo(1) % xb)
  gr_sbBodyInfo(1) % boundBox(HIGH,IAXIS) = maxval(gr_sbBodyInfo(1) % xb)
  gr_sbBodyInfo(1) % boundBox(LOW,JAXIS) = minval(gr_sbBodyInfo(1) % yb)
  gr_sbBodyInfo(1) % boundBox(HIGH,JAXIS) = maxval(gr_sbBodyInfo(1) % yb)
  gr_sbBodyInfo(1) % boundBox(LOW,KAXIS) = minval(gr_sbBodyInfo(1) % zb)
  gr_sbBodyInfo(1) % boundBox(HIGH,KAXIS) = maxval(gr_sbBodyInfo(1) % zb)

  !Creating triangles
  do i = 1, NumTriangles
     read (9, *) xb_elem(i), yb_elem(i), zb_elem(i)
     aelem(:,i) = (/xb_elem(i), yb_elem(i), zb_elem(i)/) 
  enddo

  allocate(gr_sbBodyInfo(1) % triangleCentroids(NDIM,NumTriangles))

  !Calculating centroids for each triangle
  do i = 1, NumTriangles
     gr_sbBodyInfo(1) % triangleCentroids(IAXIS,i) = &
          1./3. * (gr_sbBodyInfo(1) % xb(aelem(1,i)) + gr_sbBodyInfo(1) % xb(aelem(2,i)) + gr_sbBodyInfo(1) % xb(aelem(3,i)))

     gr_sbBodyInfo(1) % triangleCentroids(JAXIS,i) = &
          1./3. * (gr_sbBodyInfo(1) % yb(aelem(1,i)) + gr_sbBodyInfo(1) % yb(aelem(2,i)) + gr_sbBodyInfo(1) % yb(aelem(3,i)))

     gr_sbBodyInfo(1) % triangleCentroids(KAXIS,i) = &
          1./3. * (gr_sbBodyInfo(1) % zb(aelem(1,i)) + gr_sbBodyInfo(1) % zb(aelem(2,i)) + gr_sbBodyInfo(1) % zb(aelem(3,i)))
  end do

  call RuntimeParameters_get("sb_debug", gr_sbDebug)
End Subroutine gr_sbInit
