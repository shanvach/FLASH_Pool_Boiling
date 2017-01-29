!****if* source/Simulation/SimulationMain/unitTest/SolidBody/MultipleWithTriangles/MovingBodies/gr_sbInit
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
!!  * Read input file containing the coordinates of the vertices and the triangle elements of one body (ref body).
!!  * Get the min and max values of physical domain, and calculate the positions where the bodies need to be placed (in matrix form) 
!!  * Get the boundary box values
!!  * Create triangles representing each body
!!  * Calculate centroids of each triangle
!!
!! ARGUMENTS
!!

#include "constants.h"
#include "Flash.h"

Subroutine gr_sbInit()
  use Grid_data, ONLY : gr_meshMe, gr_meshComm, gr_meshNumProcs
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, &
       gr_sbPtNumX, gr_sbPtNumY, gr_sbPtNumZ, gr_sbDebug, &
       aelem, NumTriangles, NumVertices
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none
  include "Flash_mpi.h"
  

  integer :: i, j, b, k, l, ix, iy, iz, ibody, grx_sbNumBodies, gry_sbNumBodies, grz_sbNumBodies, ierr
  real, dimension(MDIM) :: body_centroid
  real, allocatable, dimension(:) :: xb_elem, yb_elem, zb_elem, xpos, ypos, zpos
  real, allocatable, dimension(:) :: ref_xb, ref_yb, ref_zb
  real, allocatable, dimension(:,:) :: ref_triangleCentroids
  integer, allocatable, dimension (:) :: proc_array
  real :: x,pd_xmin, pd_xmax, pd_ymin, pd_ymax, pd_zmin, pd_zmax
  real :: deltax, deltay, deltaz
  integer :: temp, blkID, count, myid, NProcs
  real, dimension(MDIM) :: lowOverlap, highOverlap
  real, dimension(2) :: overlapTuple
  integer,dimension(MAXBLOCKS) :: listOfBlocks
  real, dimension(2,MDIM) :: boundBox

  open (unit=9, file = 'RBC.dat', status = 'old', action = 'read')
  read (9, '(I4)') NumVertices

  read (9, '(I4)') NumTriangles

  allocate(aelem(3,NumTriangles))  !Triangles
  allocate(xb_elem(NumTriangles), yb_elem(NumTriangles), zb_elem(NumTriangles)) !Vertex numbers of trianlge

  allocate(ref_xb(NumVertices), ref_yb(NumVertices), ref_zb(NumVertices))  !Coordinate Positions of Vertices
  allocate(ref_triangleCentroids(NDIM,NumTriangles))

  call RuntimeParameters_get("sb_NumBodies", gr_sbNumBodies)
  allocate(gr_sbBodyInfo(gr_sbNumBodies))
  do b = 1, gr_sbNumBodies
     allocate(gr_sbBodyInfo(b) % xb(NumVertices), gr_sbBodyInfo(b) % yb(NumVertices), &
          gr_sbBodyInfo(b) % zb(NumVertices)) !Vertex points
     allocate(gr_sbBodyInfo(b) % triangleCentroids(NDIM,NumTriangles))
  enddo

  allocate(proc_array(0:gr_meshNumProcs-1)) !This array contains the PIds of the processors that overlap the body and 0 otherwise.

  call Grid_getListOfBlocks(LEAF, listOfBlocks, count)

  !Dividing no of bodies in the x,y,z directions
  if (NDIM == 1) then
     grx_sbNumBodies = gr_sbNumBodies
     gry_sbNumBodies = 1
     grz_sbNumBodies = 1
  endif
  if (NDIM ==2) then
     grx_sbNumBodies = gr_sbNumBodies**(1.0/2.0)
     grz_sbNumBodies = 1
     if (mod(gr_sbNumBodies,grx_sbNumBodies) == 0) then
        grx_sbNumBodies = int(gr_sbNumBodies**(1.0/2.0))
        gry_sbNumBodies = grx_sbNumBodies
     else
        grx_sbNumBodies = int(gr_sbNumBodies**(1.0/2.0))
        gry_sbNumBodies = gr_sbNumBodies - grx_sbNumBodies
     endif
  endif
  if (NDIM == 3) then
     temp = nint(gr_sbNumBodies**(1.0/3.0))
     if (temp**3 == gr_sbNumBodies) then
        grx_sbNumBodies = temp !int(gr_sbNumBodies**(1.0/3.0))
        gry_sbNumBodies = grx_sbNumBodies
        grz_sbNumBodies = grx_sbNumBodies
     else
        grx_sbNumBodies = nint(gr_sbNumBodies**(1.0/3.0))
        gry_sbNumBodies = grx_sbNumBodies
        grz_sbNumBodies= int(gr_sbNumBodies/(grx_sbNumBodies*gry_sbNumBodies))+1 
        
     endif
  endif

  allocate(xpos(gr_sbNumBodies), ypos(gr_sbNumBodies), zpos(gr_sbNumBodies)) !positions in matrix form
  xpos = 0.0
  ypos = 0.0
  zpos = 0.0
  
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

  do i = 1, NumVertices
     read (9, *) ref_xb(i), ref_yb(i), ref_zb(i) !read position coordinates of vertices (reference body)
  enddo
 
  !Creating triangles
  do i = 1, NumTriangles
     read (9, *) xb_elem(i), yb_elem(i), zb_elem(i)
     aelem(:,i) = (/xb_elem(i), yb_elem(i), zb_elem(i)/) 
  enddo
  
   !Calculating centroids for each triangle
  do i = 1, NumTriangles
     ref_triangleCentroids(IAXIS,i) = &
          1./3. * (ref_xb(aelem(1,i)) + ref_xb(aelem(2,i)) + ref_xb(aelem(3,i)))
     
     ref_triangleCentroids(JAXIS,i) = &
          1./3. * (ref_yb(aelem(1,i)) + ref_yb(aelem(2,i)) + ref_yb(aelem(3,i)))
     
     ref_triangleCentroids(KAXIS,i) = &
          1./3. * (ref_zb(aelem(1,i)) + ref_zb(aelem(2,i)) + ref_zb(aelem(3,i)))
  end do

  !Calulate centroid of the solid body = avergage of triangle centroids
  body_centroid(IAXIS) = sum(ref_triangleCentroids(IAXIS,:))/NumTriangles
  body_centroid(JAXIS) = sum(ref_triangleCentroids(JAXIS,:))/NumTriangles
  body_centroid(KAXIS) = sum(ref_triangleCentroids(KAXIS,:))/NumTriangles

  deltax = abs(pd_xmin - pd_xmax)/grx_sbNumBodies
  deltay = abs(pd_ymin - pd_ymax)/gry_sbNumBodies
  deltaz = abs(pd_zmin - pd_zmax)/grz_sbNumBodies

  !Creating many bodies
  ibody =0
  LOOP_BODY : do ix = 1, grx_sbNumBodies
     do iy = 1, gry_sbNumBodies
        do iz = 1, grz_sbNumBodies
           ibody = ibody + 1
           if (ibody > gr_sbNumBodies) exit LOOP_BODY
           xpos(ibody) = (ix-0.5)*deltax + pd_xmin
           ypos(ibody) = (iy-0.5)*deltay + pd_ymin
           zpos(ibody) = (iz-0.5)*deltaz + pd_zmin
        enddo
     enddo
  enddo LOOP_BODY

  myid = gr_meshMe 
  do b = 1, gr_sbNumBodies
     overlapTuple(1) = 0.0
     NProcs = 0
     proc_array(:) = 0

     do i = 1, NumVertices
        gr_sbBodyInfo(b)  % xb(i) = ref_xb(i) + (xpos(b) - body_centroid(IAXIS))
        gr_sbBodyInfo(b)  % yb(i) = ref_yb(i) + (ypos(b) - body_centroid(JAXIS))
        gr_sbBodyInfo(b)  % zb(i) = ref_zb(i) + (zpos(b) - body_centroid(KAXIS))
     enddo

     gr_sbBodyInfo(b) % boundBox(LOW,IAXIS) = minval(gr_sbBodyInfo(b) % xb)
     gr_sbBodyInfo(b) % boundBox(HIGH,IAXIS) = maxval(gr_sbBodyInfo(b) % xb)
     gr_sbBodyInfo(b) % boundBox(LOW,JAXIS) = minval(gr_sbBodyInfo(b) % yb)
     gr_sbBodyInfo(b) % boundBox(HIGH,JAXIS) = maxval(gr_sbBodyInfo(b) % yb)
     gr_sbBodyInfo(b) % boundBox(LOW,KAXIS) = minval(gr_sbBodyInfo(b) % zb)
     gr_sbBodyInfo(b) % boundBox(HIGH,KAXIS) = maxval(gr_sbBodyInfo(b) % zb)
           
     do i = 1, NumTriangles
        gr_sbBodyInfo(b) % triangleCentroids(IAXIS,i) = &
             1./3. * (gr_sbBodyInfo(b) % xb(aelem(1,i)) + gr_sbBodyInfo(b) % xb(aelem(2,i)) + gr_sbBodyInfo(b) % xb(aelem(3,i)))
        gr_sbBodyInfo(b) % triangleCentroids(JAXIS,i) = &
             1./3. * (gr_sbBodyInfo(b) % yb(aelem(1,i)) + gr_sbBodyInfo(b) % yb(aelem(2,i)) + gr_sbBodyInfo(b) % yb(aelem(3,i)))
        gr_sbBodyInfo(b) % triangleCentroids(KAXIS,i) = &
             1./3. * (gr_sbBodyInfo(b) % zb(aelem(1,i)) + gr_sbBodyInfo(b) % zb(aelem(2,i)) + gr_sbBodyInfo(b) % zb(aelem(3,i)))
     end do
     do i = 1, count
        blkID = listOfBlocks(i)
        call Grid_getBlkBoundBox(blkID, boundBox)
        
        if (all(&
             gr_sbBodyInfo(b) % boundBox(LOW,1:NDIM) < boundBox(HIGH,1:NDIM) .and. &
             gr_sbBodyInfo(b) % boundBox(HIGH,1:NDIM) > boundBox(LOW,1:NDIM))) then
        
           lowOverlap(1:NDIM) = max(&
                gr_sbBodyInfo(b) % boundBox(LOW,1:NDIM),&
                boundBox(LOW,1:NDIM))
           
           highOverlap(1:NDIM) = min(&
                gr_sbBodyInfo(b) % boundBox(HIGH,1:NDIM),&
                boundBox(HIGH,1:NDIM))
           
           overlapTuple(1) = overlapTuple(1) + &
                product(highOverlap(1:NDIM) - lowOverlap(1:NDIM))
          end if
     end do

!    Nprocs = 1 for each processor that overlaps. 
!    E.g. if we have four processors in total, and proc 0 and 2 overlap the body, then proc_array(for proc0) = [1 0 0 0]
!    proc_array(for proc2) = [0 0 1 0]
     if (overlapTuple(1) /= 0.0) then  
        proc_array(myid) = 1
        NProcs = NProcs + 1
     endif

     call MPI_AllReduce(NProcs, gr_sbBodyInfo(b) % nprocs, 1, &
                MPI_INT, MPI_SUM, gr_meshComm, ierr) !Total number of processors that overlap a body

     allocate(gr_sbBodyInfo(b) % communicator(gr_sbBodyInfo(b) % nprocs))

!    proc_array(in all procs) = [1 0 1 0]
     call MPI_AllReduce(MPI_IN_PLACE, proc_array, gr_meshNumProcs, &
               MPI_INT, MPI_SUM, gr_meshComm, ierr) 
     l = 1
     do k = 0, gr_meshNumProcs -1
        if (proc_array(k) > 0) then
           gr_sbBodyInfo(b) % communicator(l) = k
           l = l + 1
        end if
     end do
!     print *, "body", b, "array", gr_sbBodyInfo(b) % communicator
  enddo

  deallocate(xb_elem, yb_elem, zb_elem)
  deallocate(xpos, ypos, zpos)
  deallocate(ref_xb, ref_yb, ref_zb)
  deallocate(ref_triangleCentroids)  
  deallocate(proc_array)

  call RuntimeParameters_get("sb_debug", gr_sbDebug)
End Subroutine gr_sbInit
