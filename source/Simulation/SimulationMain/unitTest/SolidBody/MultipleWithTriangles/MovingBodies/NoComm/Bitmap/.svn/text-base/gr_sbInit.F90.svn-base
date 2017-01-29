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
!!  * Called from Grid_initDomain
!!  * Read input file containing the coordinates of the vertices and the triangle elements of one body (ref body).
!!  * Get the min and max values of physical domain, and calculate the positions where the bodies need to be placed (in matrix form) 
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
       aelem, NumAelem, NumVertices, NodesPerElem,gr_sbBodyReadProc
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none
  include "Flash_mpi.h"
  

  integer :: i, j, b, k, l, ix, iy, iz, ibody, grx_sbNumBodies, gry_sbNumBodies, grz_sbNumBodies, ierr
  real, dimension(MDIM) :: body_centroid
  real, allocatable, dimension(:) :: xb_elem, yb_elem, zb_elem, xpos, ypos, zpos
  real, allocatable, dimension(:) :: ref_xb, ref_yb, ref_zb
  real, allocatable, dimension(:,:) :: ref_triangleCentroids
  real :: x,pd_xmin, pd_xmax, pd_ymin, pd_ymax, pd_zmin, pd_zmax
  real :: deltax, deltay, deltaz
  integer :: temp, blkID, count, myid!, NProcs
  integer,dimension(MAXBLOCKS) :: listOfBlocks
  real, dimension(2,MDIM) :: boundBox

  open (unit=9, file = 'RBC.dat', status = 'old', action = 'read')
  read (9, '(I4)') NumVertices

  read (9, '(I4)') NumAelem

  NodesPerElem = 3 !Triangles

  allocate(aelem(NodesPerElem,NumAelem))  !Triangles
  allocate(xb_elem(NumAelem), yb_elem(NumAelem), zb_elem(NumAelem)) !Vertex numbers of trianlge

  allocate(ref_xb(NumVertices), ref_yb(NumVertices), ref_zb(NumVertices))  !Coordinate Positions of Vertices
  allocate(ref_triangleCentroids(NDIM,NumAelem))

  call RuntimeParameters_get("sb_NumBodies", gr_sbNumBodies)
  allocate(gr_sbBodyInfo(gr_sbNumBodies))
  allocate(gr_sbBodyReadProc(gr_sbNumBodies))

  do b = 1, gr_sbNumBodies

     body_centroid(:) = 0.0
     ref_triangleCentroids = 0.0

!     allocate(gr_sbBodyInfo(b) % xb(NumVertices), gr_sbBodyInfo(b) % yb(NumVertices), &
!          gr_sbBodyInfo(b) % zb(NumVertices)) !Vertex points
!     allocate(gr_sbBodyInfo(b) % triangleCentroids(NDIM,NumAelem))
  enddo

  call Grid_getListOfBlocks(LEAF, listOfBlocks, count)

  !Dividing no of bodies in the x,y,z directions
  if (NDIM == 1) then
     grx_sbNumBodies = gr_sbNumBodies
     gry_sbNumBodies = 1
     grz_sbNumBodies = 1
  endif
  if (NDIM ==2) then
     temp = nint(gr_sbNumBodies**(1.0/2.0))
     grz_sbNumBodies = 1
     if (temp**2 == gr_sbNumBodies) then
        grx_sbNumBodies = temp
        gry_sbNumBodies = grx_sbNumBodies
     else
        grx_sbNumBodies = nint(gr_sbNumBodies**(1.0/2.0))
        gry_sbNumBodies= int(gr_sbNumBodies/grx_sbNumBodies)+1
     endif
!     print *, grx_sbNumBodies, gry_sbNumBodies
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
     pd_ymin = 1.0
     pd_ymax = 1.0
  end if
  if (NDIM == 3) then
     call RuntimeParameters_get("sb_ptNumZ", gr_sbPtNumZ)
     call RuntimeParameters_get("zmin", pd_zmin)
     call RuntimeParameters_get("zmax", pd_zmax)
  else
     gr_sbPtNumZ = 1
     pd_zmin = 1.0
     pd_zmax = 1.0
  end if

  do i = 1, NumVertices
     read (9, *) ref_xb(i), ref_yb(i), ref_zb(i) !read position coordinates of vertices (reference body)
  enddo
 
  !Creating triangles
  do i = 1, NumAelem
     read (9, *) xb_elem(i), yb_elem(i), zb_elem(i)
     aelem(:,i) = (/xb_elem(i), yb_elem(i), zb_elem(i)/) 
  enddo
  
   !Calculating centroids for each triangle
  do i = 1, NumAelem
     ref_triangleCentroids(IAXIS,i) = &
          1./3. * (ref_xb(aelem(1,i)) + ref_xb(aelem(2,i)) + ref_xb(aelem(3,i)))
     
     if(NDIM > 1) then
        ref_triangleCentroids(JAXIS,i) = &
             1./3. * (ref_yb(aelem(1,i)) + ref_yb(aelem(2,i)) + ref_yb(aelem(3,i)))
     endif
     
     if(NDIM > 2) then
        ref_triangleCentroids(KAXIS,i) = &
             1./3. * (ref_zb(aelem(1,i)) + ref_zb(aelem(2,i)) + ref_zb(aelem(3,i)))
     endif
  end do

  !Calulate centroid of the solid body = avergage of triangle centroids
  body_centroid(IAXIS) = sum(ref_triangleCentroids(IAXIS,:))/NumAelem
  if(NDIM > 1) then
     body_centroid(JAXIS) = sum(ref_triangleCentroids(JAXIS,:))/NumAelem
  endif
  if(NDIM > 2) then
     body_centroid(KAXIS) = sum(ref_triangleCentroids(KAXIS,:))/NumAelem
  endif

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
!           print *, xpos(ibody), ypos(ibody), zpos(ibody)
        enddo
     enddo
  enddo LOOP_BODY

 ! This read is hacked for now - all bodies read by MASTER_PE - MV
  gr_sbBodyReadProc(1:gr_sbNumBodies) = MASTER_PE !Ques: Where is this routine? More about this routine.

  do b = 1, gr_sbNumBodies

     gr_sbBodyInfo(b)%bodyMaster = gr_sbBodyReadProc(b)

     gr_sbBodyInfo(b)%NumVertices = NumVertices
     gr_sbBodyInfo(b)%NumAelem = NumAelem

     if (gr_meshMe .eq. gr_sbBodyInfo(b)%bodyMaster) then 

     ! Then allocate the other fields in gr_sbBodyInfo
        NumVertices  = gr_sbBodyInfo(b)%NumVertices
        NumAelem = gr_sbBodyInfo(b)%NumAelem

        allocate(gr_sbBodyInfo(b) % xbus(NumVertices), gr_sbBodyInfo(b) % ybus(NumVertices))
        allocate(gr_sbBodyInfo(b) % xb(NumVertices), gr_sbBodyInfo(b) % yb(NumVertices)) !Vertex points
        allocate(gr_sbBodyInfo(b) % ubd(NumVertices), gr_sbBodyInfo(b) % vbd(NumVertices)) !Vertex velocities
        allocate(gr_sbBodyInfo(b) % ubdd(NumVertices), gr_sbBodyInfo(b) % vbdd(NumVertices)) !Vertex accelrations
        allocate(gr_sbBodyInfo(b) % AAREANODE(NumVertices), gr_sbBodyInfo(b) % AANGNODE(NumVertices))
        allocate(gr_sbBodyInfo(b) % dutdn(NumVertices), gr_sbBodyInfo(b) % dvtdn(NumVertices))
        allocate(gr_sbBodyInfo(b) % nxL(NumVertices), gr_sbBodyInfo(b) % nyL(NumVertices))   !  nzL for 3D
        allocate(gr_sbBodyInfo(b) % AAREAELEM(NumAelem))
        allocate(gr_sbBodyinfo(b) % AELEMNODE(NumVertices,10))    ! Array of elements connected to each node 
        allocate(gr_sbBodyInfo(b) % AELEM(NodesPerElem + CONSTANT_ONE,NumAelem))        ! Array of connectivity for each body
        allocate(gr_sbBodyInfo(b) % sbPtNumXi(NumAelem),gr_sbBodyInfo(b) % sbPtNumEta(NumAelem))
        allocate(gr_sbBodyInfo(b) % sbPtNumElem(NumAelem))
        allocate(gr_sbBodyInfo(b) % AELTYPE(NumAelem))

#if NDIM == 2
        allocate(gr_sbBodyInfo(b) % sb(NumVertices))
#endif

#if NDIM == 3
        allocate(gr_sbBodyInfo(b) % zbus(NumVertices))
        allocate(gr_sbBodyInfo(b) % zb(NumVertices))
        allocate(gr_sbBodyInfo(b) % wbd(NumVertices))  ! vertex velocity
        allocate(gr_sbBodyInfo(b) % wbdd(NumVertices)) ! vertex accelaration
        allocate(gr_sbBodyInfo(b) % dwtdn(NumVertices))

        allocate(gr_sbBodyInfo(b) % nzL(NumAelem))
#endif


        allocate(gr_sbBodyInfo(b) % triangleCentroids(NDIM,NumAelem))


     !Fill Bodies xb, yb, zb
        do i = 1, NumVertices
           gr_sbBodyInfo(b)  % xb(i) = ref_xb(i)*deltax + (xpos(b) - body_centroid(IAXIS))
           if (NDIM >= 2) then
              gr_sbBodyInfo(b)  % yb(i) = ref_yb(i)*deltay + (ypos(b) - body_centroid(JAXIS))
           else
              gr_sbBodyInfo(b)  % yb(i) = NONEXISTENT
           endif
           if(NDIM == 3) then
              gr_sbBodyInfo(b)  % zb(i) = ref_zb(i)*deltaz + (zpos(b) - body_centroid(KAXIS))
           else
              gr_sbBodyInfo(b)  % zb(i) = NONEXISTENT
           endif
           if ((gr_sbBodyInfo(b)  % xb(i) < pd_xmin) .or. (gr_sbBodyInfo(b)  % xb(i) > pd_xmax)) then
              print *, "body", b, "x outside domain", gr_sbBodyInfo(b)  % xb(i)
           endif
           if(NDIM > 1) then
              if ((gr_sbBodyInfo(b)  % yb(i) < pd_ymin) .or. (gr_sbBodyInfo(b)  % yb(i) > pd_ymax)) then
                 print *, "body", b, "y outside domain", gr_sbBodyInfo(b)  % yb(i)
              endif
           endif
           if(NDIM > 2) then
              if ((gr_sbBodyInfo(b)  % zb(i) < pd_zmin) .or. (gr_sbBodyInfo(b)  % zb(i) > pd_zmax)) then
                 print *, "body", b, "z outside domain", gr_sbBodyInfo(b)  % zb(i)
              endif
           endif
        enddo

#if NDIM ==2
        write(*,*) 'b=',b,NodesPerElem
        gr_sbBodyInfo(b) % AELEM(CONSTANT_ONE,1:NumAelem) = CONSTANT_TWO ! Number of vertices on a segment 
        gr_sbBodyInfo(b) % AELEM(CONSTANT_TWO:NDIM+1,1:NumAelem) = aelem(1:NDIM,1:NumAelem)
#elif NDIM == 3
        gr_sbBodyInfo(b) % zbus(1:NumVertices) = ref_zb(1:NumVertices)
        gr_sbBodyInfo(b) % AELEM(CONSTANT_ONE,1:NumAelem) = NDIM ! Number of vertices on a triangle 
        gr_sbBodyInfo(b) % AELEM(CONSTANT_TWO:NDIM+1,1:NumAelem) = aelem(1:NDIM,1:NumAelem)
#endif

     endif
  enddo

  deallocate(xb_elem, yb_elem, zb_elem)
  deallocate(xpos, ypos, zpos)
  deallocate(ref_xb, ref_yb, ref_zb)
  deallocate(ref_triangleCentroids)  
  
  call RuntimeParameters_get("sb_debug", gr_sbDebug)
End Subroutine gr_sbInit
