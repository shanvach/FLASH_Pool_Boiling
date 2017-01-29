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
       gr_sbDebug, gr_sbBodyReadProc, &
       aelem, NumAelem, NumVertices, NodesPerElem
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none
  include "Flash_mpi.h"
  

  integer :: i, j, b, k, l, ix, iy, iz, ibody, grx_sbNumBodies, gry_sbNumBodies, grz_sbNumBodies, ierr
  real, dimension(MDIM) :: body_centroid
  real, allocatable, dimension(:) :: xb_elem, yb_elem, zb_elem, xpos, ypos, zpos
  real, allocatable, dimension(:) :: ref_xb, ref_yb, ref_zb
  real, allocatable, dimension(:,:) :: ref_triangleCentroids
!  integer, allocatable, dimension (:) :: proc_array
  real :: x,pd_xmin, pd_xmax, pd_ymin, pd_ymax, pd_zmin, pd_zmax
  real :: deltax, deltay, deltaz
  integer :: temp, blkID, count, myid!, NProcs
!  real, dimension(MDIM) :: lowOverlap, highOverlap
!  real, dimension(2) :: overlapTuple
  integer,dimension(MAXBLOCKS) :: listOfBlocks
  real, dimension(2,MDIM) :: boundBox


  call RuntimeParameters_get("sb_NumBodies", gr_sbNumBodies)
  allocate(gr_sbBodyInfo(gr_sbNumBodies))
  allocate(gr_sbBodyReadProc(gr_sbNumBodies))

!  ! Read Files: 
!  open (unit=9, file = 'CYL.dat', status = 'old', action = 'read')
!  read (9, '(I8)') NumVertices
!  read (9, '(2I8)') NumAelem,NodesPerElem

  ! Read Files: 
  open (unit=9, file = 'sphere_coarse_altered.surf', status = 'old', action = 'read')
  read (9, '(I8)') NumVertices
  read (9,*) NumAelem,NodesPerElem

  ! This read is hacked for now - all bodies read by MASTER_PE - MV
  gr_sbBodyReadProc(1:gr_sbNumBodies) = MASTER_PE


  ! First Set the number of vertices and triangles:
  do b = 1, gr_sbNumBodies

     gr_sbBodyInfo(b)%bodyMaster = gr_sbBodyReadProc(b)

     gr_sbBodyInfo(b)%NumVertices  = NumVertices
     gr_sbBodyInfo(b)%NumAelem     = NumAelem


     if (gr_meshMe .eq. gr_sbBodyInfo(b)%bodyMaster) then

     ! Then allocate the other fields in gr_sbBodyInfo
     NumVertices  = gr_sbBodyInfo(b)%NumVertices 
     NumAelem = gr_sbBodyInfo(b)%NumAelem 

     allocate(aelem(MDIM,NumAelem))                                           !Aerodynamic Elements
     allocate(xb_elem(NumAelem), yb_elem(NumAelem), zb_elem(NumAelem))        !Vertex numbers in Elements
     allocate(ref_xb(NumVertices), ref_yb(NumVertices), ref_zb(NumVertices))  !Coordinate Positions of Vertices



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


     ! After allocation Fill bodies structures xbus, ybus, zbus and AELEM
     ! ----------------------------------------------------------------------------------------------------
     ! For Now read positions from same file:
     ! Vertices:
#if NDIM == 2
     do i = 1, NumVertices
        read (9, *) ref_xb(i), ref_yb(i) !read position coordinates of vertices (reference body)
     enddo
     !Creating triangles
     do i = 1, NumAelem
        read (9, *) xb_elem(i), yb_elem(i)
        aelem(1:NDIM,i) = (/xb_elem(i), yb_elem(i) /) 
     enddo
#elif NDIM == 3

     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     do i = 1, NumVertices
        read (9, *) ref_xb(i), ref_yb(i), ref_zb(i) !read position coordinates of vertices (reference body)
     enddo
     !Creating triangles
     do i = 1, NumAelem
        read (9, *) xb_elem(i), yb_elem(i), zb_elem(i)
        aelem(1:NDIM,i) = (/xb_elem(i), yb_elem(i), zb_elem(i)/) 
     enddo
     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#endif
     ! ----------------------------------------------------------------------------------------------------
     write(*,*) 'b=',b,' ,NumBods=',gr_sbNumBodies

     ! Load into gr_sbBodyInfo
     gr_sbBodyInfo(b) % xbus(1:NumVertices) = ref_xb(1:NumVertices)
     gr_sbBodyInfo(b) % ybus(1:NumVertices) = ref_yb(1:NumVertices)
#if NDIM ==2
     write(*,*) 'b=',b,NodesPerElem
     gr_sbBodyInfo(b) % AELEM(CONSTANT_ONE,1:NumAelem) = CONSTANT_TWO ! Number of vertices on a segment 
     gr_sbBodyInfo(b) % AELEM(CONSTANT_TWO:NDIM+1,1:NumAelem) = aelem(1:NDIM,1:NumAelem)
#elif NDIM == 3
     write(*,*) 'b=',b,NodesPerElem
     gr_sbBodyInfo(b) % zbus(1:NumVertices) = ref_zb(1:NumVertices)
     gr_sbBodyInfo(b) % AELEM(CONSTANT_ONE,1:NumAelem) = NDIM ! Number of vertices on a triangle 
     gr_sbBodyInfo(b) % AELEM(CONSTANT_TWO:NDIM+1,1:NumAelem) = aelem(1:NDIM,1:NumAelem)     
#endif
     

     deallocate(xb_elem, yb_elem, zb_elem)
     deallocate(ref_xb, ref_yb, ref_zb)
     deallocate(aelem)

  endif

  enddo


  call RuntimeParameters_get("sb_debug", gr_sbDebug)
End Subroutine gr_sbInit
