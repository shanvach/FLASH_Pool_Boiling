!!****if* source/Grid/GridStructures/withTriangles/pointTopoint/bitmap/gr_sbData
!!
!! NAME
!!  gr_sbData
!!
!! SYNOPSIS
!!
!!  use gr_sbData
!!
!!
!!***

#include "constants.h"
#include "Flash.h"

Module gr_sbData
implicit none

type solid_body
   real, pointer, dimension(:,:) :: particles
   real, pointer, dimension(:,:) :: triangleCentroids
   integer :: myPE
   integer :: bodyMaster
   integer :: sendProcs !count of processors that need to recieve particles
   real, dimension(2,MDIM) :: boundBox

   integer :: model, memflag, rigflexflag
   integer :: NumAelem,NumVertices,totalPart
   real :: xo,yo,zo
   real, allocatable, dimension(:) :: xbus, ybus, zbus
   real, allocatable, dimension(:) :: xb, yb, zb, sb !Vertex coordinates
   real, allocatable, dimension(:) :: ubd, vbd, wbd
   real, allocatable, dimension(:) :: ubdd, vbdd, wbdd
   real, allocatable, dimension(:) :: nxL, nyL, nzL
   real, allocatable, dimension(:) :: dutdn, dvtdn, dwtdn
   integer, allocatable, dimension(:,:) :: AELEMNODE      ! Array of aerodynamic elements conected to the aerodynamic nodes.   
   real, allocatable, dimension(:)      :: AAREANODE,AANGNODE ! Array of node associated areas
   real, allocatable, dimension(:)      :: AAREAELEM
   integer, allocatable, dimension(:,:) :: AELEM
   integer, allocatable, dimension(:)   :: sbPtNumXi, sbPtNumEta
   integer, allocatable, dimension(:)   :: sbPtNumElem 
   integer, allocatable, dimension(:)   :: AELTYPE

   integer, pointer, dimension(:,:) :: particlesPerProc !contains particle count of each proc. This information contained in master
   integer :: sbIsFixed   = CONSTANT_ZERO  ! Body, by default is not fixed.

   integer, allocatable, dimension(:,:,:,:) :: ielem
   real,    allocatable, dimension(:,:,:,:) :: phile 

end type solid_body

integer, allocatable, dimension(:) :: gr_sbParticleCount !Particle count stored in appropriate processors

integer, allocatable, dimension(:) :: gr_sbBodyReadProc

real, allocatable, dimension(:,:) :: aelem !Aerodynamic grid elements

type(solid_body), save, dimension(:), pointer :: gr_sbBodyInfo
integer, save :: gr_sbNumBodies    ! total number of bodies in the computation

integer, save :: NodesPerElem
integer, save :: NumAelem
integer, save :: NumVertices

integer, save :: totalPart, sumPart !sumPart: sum of all particles from all bodies a proc has

!Number of particles to generate within each triangle.
integer, save :: gr_sbPtNumXi, gr_sbPtNumX, gr_sbPtNumY, gr_sbPtNumZ
integer, save :: gr_sbPtNumEta

logical, save :: gr_sbDebug

integer, save, allocatable, dimension(:,:,:) :: gr_sbIJKtoProc ! Variable only used in uniform grids.

integer, save :: gr_sbFirstCall = CONSTANT_ONE   ! First call is initialized true.

integer, save :: gr_sbStencil

End Module gr_sbData
