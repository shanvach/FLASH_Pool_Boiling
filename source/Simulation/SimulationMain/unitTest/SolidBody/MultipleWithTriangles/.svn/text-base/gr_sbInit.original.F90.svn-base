!Called from Grid_init
!Allocate and populate the data structure that holds all information
!about the single solid body.

#include "constants.h"
#include "Flash.h"

Subroutine gr_sbInit()
  use Grid_data, ONLY : gr_meshMe
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, &
       gr_sbPtNumX, gr_sbPtNumY, gr_sbPtNumZ, gr_sbDebug, &
       xb, yb, zb, aelem
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none
  integer, parameter :: NumTriangles = 4
  integer, parameter :: NumVertices = 6
  integer :: i

  allocate(xb(NumVertices), yb(NumVertices), zb(NumVertices)) !Vertex points
  allocate(aelem(3,NumTriangles))  !Triangles

  gr_sbNumBodies = 1
  allocate(gr_sbBodyInfo(gr_sbNumBodies))

  gr_sbBodyInfo(1) % boundBox(:,:) = 0.0

  call RuntimeParameters_get("sb_xmin", gr_sbBodyInfo(1) % boundBox(LOW,IAXIS))
  call RuntimeParameters_get("sb_xmax", gr_sbBodyInfo(1) % boundBox(HIGH,IAXIS))
  call RuntimeParameters_get("sb_ptNumX", gr_sbPtNumX)

  if (NDIM >= 2) then
     call RuntimeParameters_get("sb_ymin", gr_sbBodyInfo(1) % boundBox(LOW,JAXIS))
     call RuntimeParameters_get("sb_ymax", gr_sbBodyInfo(1) % boundBox(HIGH,JAXIS))
     call RuntimeParameters_get("sb_ptNumY", gr_sbPtNumY)
  else
     gr_sbPtNumY = 1
  end if

  if (NDIM == 3) then
     call RuntimeParameters_get("sb_zmin", gr_sbBodyInfo(1) % boundBox(LOW,KAXIS))
     call RuntimeParameters_get("sb_zmax", gr_sbBodyInfo(1) % boundBox(HIGH,KAXIS))
     call RuntimeParameters_get("sb_ptNumZ", gr_sbPtNumZ)
  else
     gr_sbPtNumZ = 1
  end if

  xb(1) = gr_sbBodyInfo(1) % boundBox(LOW,IAXIS)
  xb(2) = gr_sbBodyInfo(1) % boundBox(HIGH,IAXIS)
  xb(3) = gr_sbBodyInfo(1) % boundBox(HIGH,IAXIS)
  xb(4) = gr_sbBodyInfo(1) % boundBox(LOW,IAXIS)
  xb(5) = gr_sbBodyInfo(1) % boundBox(HIGH,IAXIS)
  xb(6) = gr_sbBodyInfo(1) % boundBox(HIGH,IAXIS)

  yb(1) = gr_sbBodyInfo(1) % boundBox(LOW,JAXIS)
  yb(2) = gr_sbBodyInfo(1) % boundBox(LOW,JAXIS)
  yb(3) = gr_sbBodyInfo(1) % boundBox(LOW,JAXIS)
  yb(4) = gr_sbBodyInfo(1) % boundBox(LOW,JAXIS)
  yb(5) = gr_sbBodyInfo(1) % boundBox(HIGH,JAXIS)
  yb(6) = gr_sbBodyInfo(1) % boundBox(HIGH,JAXIS)

  zb(1) = gr_sbBodyInfo(1) % boundBox(LOW,KAXIS)
  zb(2) = gr_sbBodyInfo(1) % boundBox(LOW,KAXIS)
  zb(3) = gr_sbBodyInfo(1) % boundBox(HIGH,KAXIS)
  zb(4) = gr_sbBodyInfo(1) % boundBox(HIGH,KAXIS)
  zb(5) = gr_sbBodyInfo(1) % boundBox(LOW,KAXIS)
  zb(6) = gr_sbBodyInfo(1) % boundBox(HIGH,KAXIS)

  !Creating triangles
  aelem(:,1) = (/1,3,4/)
  aelem(:,2) = (/1,2,3/)
  aelem(:,3) = (/2,5,3/)
  aelem(:,4) = (/3,5,6/)

  allocate(gr_sbBodyInfo(1) % triangleCentroids(NDIM,NumTriangles))

  !Calculating centroids for each triangle
  do i = 1, NumTriangles
     gr_sbBodyInfo(1) % triangleCentroids(IAXIS,i) = &
          1./3. * (xb(aelem(1,i)) + xb(aelem(2,i)) + xb(aelem(3,i)))

     gr_sbBodyInfo(1) % triangleCentroids(JAXIS,i) = &
          1./3. * (yb(aelem(1,i)) + yb(aelem(2,i)) + yb(aelem(3,i)))

     gr_sbBodyInfo(1) % triangleCentroids(KAXIS,i) = &
          1./3. * (zb(aelem(1,i)) + zb(aelem(2,i)) + zb(aelem(3,i)))

     if (gr_meshMe == 0) then
        write(6,'(a,i4,a,3e12.2)') " Triangle ", i, " has centroid :", &
             gr_sbBodyInfo(1) % triangleCentroids(:,i)
     end if
  end do

  call RuntimeParameters_get("sb_debug", gr_sbDebug)
End Subroutine gr_sbInit
