!!****if* source/flashUtilities/contourSurface/ut_contourSurfaceAreaBlock
!!
!! NAME
!!
!!  sim_contour_surface_area
!!
!!
!! SYNOPSIS
!!
!!  contour_surface_area(ivar,level,area)
!!
!!  contour_surface_area(integer,real,real)
!!
!!
!! DESCRIPTION
!!
!!  Extracts contour surface of a given variable using marching
!!  cubes algorithm and computs area of this surface for all blocks
!!  on this processor.
!!
!!  note this is just a wrapper routine for a marching cubes routine
!!  the vertex data and box dimensions are calculated and passed to
!!  a subroutine that actually does the marching cubes.
!!
!!  this function should not be called in a 1d simulation
!!
!!
!! ARGUMENTS
!!
!!  ivar        index of variable to
!!
!!  level       contour level to process
!!
!!  area        area of the contour surface
!!
!!***

subroutine ut_contourSurfaceAreaBlock(nlevels,isolevels,ctrData, blkLimits, blkIndex, areas)

use Grid_interface, ONLY : Grid_getBlkBoundBox
use Driver_interface, ONLY : Driver_abortFlash

  implicit none
#include "Flash.h"
#include "constants.h"

  integer,                     intent(IN)  :: nlevels
  real,dimension(nlevels),     intent(IN)  :: isolevels
  real,dimension(:,:,:),       intent(IN)  :: ctrData
  integer,dimension(HIGH,MDIM),intent(IN)  :: blkLimits
  integer,                     intent(IN)  :: blkIndex
  real,dimension(nlevels),     intent(OUT) :: areas

  real, allocatable, dimension(:,:,:) :: vertData

  integer, dimension(3) :: vdataSize

  integer :: i,j,k,ii,jj,kk, istat
  real :: xmin,ymin,zmin,xmax,ymax,zmax
  real :: vdatamin, vdatamax
  real, dimension(2,MDIM) :: boundBox

  call Grid_getBlkBoundBox(blkIndex, boundBox)

  xmin = boundBox(1,IAXIS)
  ymin = boundBox(1,JAXIS)
  zmin = boundBox(1,KAXIS)
  xmax = boundBox(2,IAXIS)
  ymax = boundBox(2,JAXIS)
  zmax = boundBox(2,KAXIS)
#if NDIM == 2
  zmax = zmin+xmax-xmin
#endif

  ! allocate space for vertex data
  ! always keep two cells in z direction for faking
  vdataSize(1) = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+2
  vdataSize(2) = blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+2
  vdataSize(3) = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+2
  allocate(vertData(vdataSize(1), vdataSize(2), vdataSize(3)),STAT=istat)
   if (istat /= 0) call Driver_abortFlash("Cannot allocate vertdata in sim_contour_surface_area")

  vdatamax = -HUGE(1.0)
  vdatamin = HUGE(1.0)

  ! calculate vertex data by averaging neighboring cells
#if NDIM == 2
  kk = blkLimits(LOW,KAXIS)
  do j = 1, vdataSize(2)
     jj = blkLimits(LOW,JAXIS)-2 + j
     do i = 1, vdataSize(1)
        ii = blkLimits(LOW,IAXIS)-2 + i
        vertData(i,j,1) = 0.25*( ctrData(ii,jj,kk)   + ctrData(ii+1,jj,kk) + &
                                 ctrData(ii,jj+1,kk) + ctrData(ii+1,jj+1,kk) )
        vertData(i,j,2) = vertdata(i,j,1)
        vdatamax = max(vdatamax,vertdata(i,j,1))
        vdatamin = min(vdatamin,vertdata(i,j,1))
     end do
  end do
#endif

#if NDIM == 3
  do k = 1, vdataSize(3)
     kk = blkLimits(LOW,KAXIS)-2 + k
     do j = 1, vdataSize(2)
        jj = blkLimits(LOW,JAXIS)-2 + j
        do i = 1, vdataSize(1)
           ii = blkLimits(LOW,IAXIS)-2 + i
           vertData(i,j,k) = 0.125*( ctrData(ii,jj,kk)     + ctrData(ii+1,jj,kk) + &
                                     ctrData(ii,jj+1,kk)   + ctrData(ii+1,jj+1,kk) + &
                                     ctrData(ii,jj,kk+1)   + ctrData(ii+1,jj,kk+1) + &
                                     ctrData(ii,jj+1,kk+1) + ctrData(ii+1,jj+1,kk+1) )
           vdatamax = max(vdatamax,vertData(i,j,k))
           vdatamin = min(vdatamin,vertData(i,j,k))
        end do
     end do
  end do
#endif

  do i = 1, nlevels
     if( vdatamin <= isolevels(i) .AND. isolevels(i) <= vdatamax ) then
        call iso_surface(vertData,vdatasize(1),vdataSize(2),vdataSize(3), &
                         xmin,ymin,zmin,xmax,ymax,zmax,isolevels(i),areas(i))
#if NDIM == 2
        areas(i) = areas(i)/(zmax-zmin)
#endif
     else
        areas(i) = 0.0
     end if
  enddo

  deallocate(vertData)

end subroutine ut_contourSurfaceAreaBlock


