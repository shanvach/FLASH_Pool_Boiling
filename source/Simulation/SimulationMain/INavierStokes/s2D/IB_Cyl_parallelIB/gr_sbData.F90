!!****if* source/Simulation/SimulationMain/INavierStokes/2D/IB_Cyl_parallelIB/gr_sbData
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

Module gr_sbData
implicit none

type solid_body
   real, pointer, dimension(:,:) :: particles
   real, pointer, dimension(:,:) :: triangleCentroids
   integer :: myPE
   integer :: bodyMaster
   integer :: sendProcs
   real, dimension(2,MDIM) :: boundBox
   integer :: model, memflag, rigflexflag
   integer :: NumTriangles,NumVertices,totalPart
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
   integer, allocatable, dimension(:)   :: sbPtNumX, sbPtNumY, sbPtNumZ 
   integer, pointer, dimension(:,:) :: particlesPerProc !contains particle count of each proc. This information contained in master
end type solid_body

integer, allocatable, dimension(:) :: gr_sbParticleCount !Particle count stored in appropriate processors

real, allocatable, dimension(:,:) :: aelem !Triangle elements

type(solid_body), save, dimension(:), pointer :: gr_sbBodyInfo
integer, save :: gr_sbNumBodies    ! total number of bodies in the computation

integer, save :: NodesPerElem
integer, save :: NumTriangles
integer, save :: NumVertices

integer, save :: totalPart !total no of particles for each body

!Number of particles to generate within each triangle.
integer, save :: gr_sbPtNumX
integer, save :: gr_sbPtNumY
integer, save :: gr_sbPtNumZ

logical, save :: gr_sbDebug

End Module gr_sbData
