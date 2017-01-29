!****if* source/Simulation/SimulationMain/unitTest/SolidBody/MultipleWithTriangles/gr_sbInit
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
!!  * Read input file containing the coordinates of the vertices and the triangle elements of one body.
!!  * Have many bodies. The position of the multiple bodies is randomly generated. The position is accepted 
!!     if no overlap with bodies already placed (Acceptance criterion). The body to be placed is compared with 
!!     all bodies, n-1 comaprisons.
!!  * Get the boundary box values
!!  * Create triangles representing each body
!!  * Calculate centroids of each triangle
!!
!! ARGUMENTS
!!

#include "constants.h"
#include "Flash.h"

Subroutine gr_sbInit()
  use Grid_data, ONLY : gr_meshComm
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, &
       gr_sbPtNumX, gr_sbPtNumY, gr_sbPtNumZ, gr_sbDebug, &
       aelem, NumTriangles, NumVertices
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none
  include "Flash_mpi.h"
  

  integer :: ierr, id, i, j, b, loopcount
  real, allocatable, dimension(:) :: xb_elem, yb_elem, zb_elem, temp_xb, temp_yb, temp_zb
  real :: x,shift, pd_xmin, pd_xmax, pd_ymin, pd_ymax, pd_zmin, pd_zmax

  open (unit=9, file = 'RBC.dat', status = 'old', action = 'read')
  read (9, '(I4)') NumVertices

  read (9, '(I4)') NumTriangles

  allocate(aelem(3,NumTriangles))  !Triangles
  allocate(xb_elem(NumTriangles), yb_elem(NumTriangles), zb_elem(NumTriangles)) !Vertex numbers of trianlge

  call RuntimeParameters_get("sb_NumBodies", gr_sbNumBodies)
  allocate(gr_sbBodyInfo(gr_sbNumBodies))

  do b = 1, gr_sbNumBodies
     allocate(gr_sbBodyInfo(b) % xb(NumVertices), gr_sbBodyInfo(b) % yb(NumVertices), &
          gr_sbBodyInfo(b) % zb(NumVertices)) !Vertex points
     allocate(gr_sbBodyInfo(b) % triangleCentroids(NDIM,NumTriangles))
  enddo

  allocate(temp_xb(NumVertices), temp_yb(NumVertices), temp_zb(NumVertices))

  temp_xb = 0
  temp_yb = 0
  temp_zb = 0
  
  !No of particles to be created in each triangle and the physical domain
  call RuntimeParameters_get("sb_ptNumX", gr_sbPtNumX)
  call RuntimeParameters_get("xmin", pd_xmin)
  call RuntimeParameters_get("xmax", pd_xmax)
  if (NDIM >= 2) then
     call RuntimeParameters_get("sb_ptNumY", gr_sbPtNumY)
     call RuntimeParameters_get("ymin", pd_ymin)
     call RuntimeParameters_get("ymax", pd_ymax)
  else
     gr_sbPtNumY = 1
  end if
  if (NDIM == 3) then
     call RuntimeParameters_get("sb_ptNumZ", gr_sbPtNumZ)
     call RuntimeParameters_get("zmin", pd_zmin)
     call RuntimeParameters_get("zmax", pd_zmax)
  else
     gr_sbPtNumZ = 1
  end if

  gr_sbBodyInfo(1) % boundBox(:,:) = 0.0
  
  do i = 1, NumVertices
     read (9, *) gr_sbBodyInfo(1) % xb(i), gr_sbBodyInfo(1) % yb(i), gr_sbBodyInfo(1) % zb(i) !read position coordinates of vertices
  enddo
 
  !Boundary box values of the body
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
  
   !Calculating centroids for each triangle
  do i = 1, NumTriangles
     gr_sbBodyInfo(1) % triangleCentroids(IAXIS,i) = &
          1./3. * (gr_sbBodyInfo(1) % xb(aelem(1,i)) + gr_sbBodyInfo(1) % xb(aelem(2,i)) + gr_sbBodyInfo(1) % xb(aelem(3,i)))
     
     gr_sbBodyInfo(1) % triangleCentroids(JAXIS,i) = &
          1./3. * (gr_sbBodyInfo(1) % yb(aelem(1,i)) + gr_sbBodyInfo(1) % yb(aelem(2,i)) + gr_sbBodyInfo(1) % yb(aelem(3,i)))
     
     gr_sbBodyInfo(1) % triangleCentroids(KAXIS,i) = &
          1./3. * (gr_sbBodyInfo(1) % zb(aelem(1,i)) + gr_sbBodyInfo(1) % zb(aelem(2,i)) + gr_sbBodyInfo(1) % zb(aelem(3,i)))
  end do
  
  !Creating many bodies, random positions (Monte Carlo Simulation)
  if (gr_sbNumBodies > 1) then
     call random_seed()
     do b = 2, gr_sbNumBodies
        temp_xb = 0.0
        temp_yb = 0.0
        temp_zb = 0.0

     !   gr_sbBodyInfo(b) % boundBox(:,:) = 0.0
        loopcount = 0
    12    call MPI_Comm_rank(gr_meshComm, id, ierr)
        !Random placement of bodies
        if (id == 0) then
           gr_sbBodyInfo(b) % boundBox(:,:) = 0.0
           call random_number(x)
           shift = (pd_xmin - minval(gr_sbBodyInfo(1) % xb))+((pd_xmax - maxval(gr_sbBodyInfo(1) % xb)) - (pd_xmin - minval(gr_sbBodyInfo(1) % xb)))*x
        endif

        call MPI_BCAST(shift, 1, FLASH_REAL, 0, gr_meshComm, ierr)

        do i = 1, NumVertices
           temp_xb(i) = gr_sbBodyInfo(1)  % xb(i) + shift
           temp_yb(i) = gr_sbBodyInfo(1)  % yb(i) + shift
           temp_zb(i) = gr_sbBodyInfo(1)  % zb(i) + shift
        enddo
        
        gr_sbBodyInfo(b) % boundBox(LOW,IAXIS) = minval(temp_xb)
        gr_sbBodyInfo(b) % boundBox(HIGH,IAXIS) = maxval(temp_xb)
        gr_sbBodyInfo(b) % boundBox(LOW,JAXIS) = minval(temp_yb)
        gr_sbBodyInfo(b) % boundBox(HIGH,JAXIS) = maxval(temp_yb)
        gr_sbBodyInfo(b) % boundBox(LOW,KAXIS) = minval(temp_zb)
        gr_sbBodyInfo(b) % boundBox(HIGH,KAXIS) = maxval(temp_zb)

        !Acceptance criteria
        do j = 1, b-1
           if (all(&
             gr_sbBodyInfo(b) % boundBox(LOW,1:NDIM) < gr_sbBodyInfo(j) % boundBox(HIGH,1:NDIM) .and. &
             gr_sbBodyInfo(b) % boundBox(HIGH,1:NDIM) > gr_sbBodyInfo(j) % boundBox(LOW,1:NDIM))) then
              loopcount = loopcount + 1
              if (loopcount < 100) then  !to check the body not to go in an infinite loop, will change when comparing only with neighboring bodies.
                 go to 12
              endif
           else
              if (j == b - 1) then !The new body has been compared with all the bodies
                 gr_sbBodyInfo(b)  % xb = temp_xb
                 gr_sbBodyInfo(b)  % yb = temp_yb
                 gr_sbBodyInfo(b)  % zb = temp_zb
                 do i = 1, NumTriangles
                    gr_sbBodyInfo(b) % triangleCentroids(IAXIS,:) = &
                         1./3. * (gr_sbBodyInfo(b)  % xb(aelem(1,i)) + gr_sbBodyInfo(b)  % xb(aelem(2,i)) + gr_sbBodyInfo(b)  % xb(aelem(3,i)))
                    gr_sbBodyInfo(b) % triangleCentroids(JAXIS,:) = &
                         1./3. * (gr_sbBodyInfo(b)  % yb(aelem(1,i)) + gr_sbBodyInfo(b)  % yb(aelem(2,i)) + gr_sbBodyInfo(b)  % yb(aelem(3,i)))
                    gr_sbBodyInfo(b) % triangleCentroids(KAXIS,:) = &
                         1./3. * (gr_sbBodyInfo(b)  % zb(aelem(1,i)) + gr_sbBodyInfo(b)  % zb(aelem(2,i)) + gr_sbBodyInfo(b)  % zb(aelem(3,i)))
                 end do
                 gr_sbBodyInfo(b) % boundBox(LOW,IAXIS) = minval(gr_sbBodyInfo(b) % xb)
                 gr_sbBodyInfo(b) % boundBox(HIGH,IAXIS) = maxval(gr_sbBodyInfo(b) % xb)
                 gr_sbBodyInfo(b) % boundBox(LOW,JAXIS) = minval(gr_sbBodyInfo(b)% yb)
                 gr_sbBodyInfo(b) % boundBox(HIGH,JAXIS) = maxval(gr_sbBodyInfo(b) % yb)
                 gr_sbBodyInfo(b) % boundBox(LOW,KAXIS) = minval(gr_sbBodyInfo(b) % zb)
                 gr_sbBodyInfo(b) % boundBox(HIGH,KAXIS) = maxval(gr_sbBodyInfo(b) % zb)
                 exit
              endif
           endif
        enddo
     enddo
  endif

  deallocate(xb_elem, yb_elem, zb_elem)
  deallocate(temp_xb, temp_yb, temp_zb)
  
  call RuntimeParameters_get("sb_debug", gr_sbDebug)
End Subroutine gr_sbInit
