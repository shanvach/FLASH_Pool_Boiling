!!****if* source/Grid/GridMain/UG/Grid_colocateFaceData
!!
!! NAME
!!
!!  Grid_colocateFaceData
!!
!!
!! SYNOPSIS
!!
!!   call Grid_colocateFaceData(integer(IN) :: faceVar,
!!                              integer(IN) :: axis,
!!                              integer(IN) :: blockid,
!!                              real(OUT)   :: cellData(:,:,:))
!!
!! DESCRIPTION
!!
!! This function colocates the face centered grid variables by interpolating 
!! the face data to the cell centers.
!!
!! Currently only used for colocated velocity data in plot files for hdf5 IO
!!           and only for the uniform grid, and derived, grids
!!
!! ARGUMENTS
!!
!!  faceVar - Specifies the face variable from which to interpolate the data
!!
!!  axis - Specifies axis to which the face data belongs
!!
!!  blockid - Specifies the block for which interpolate the face data
!!
!!  cellData - The array holding the data returning the interpolated face data
!!
!!***
 
subroutine Grid_colocateFaceData(faceVar, axis, blockid, cellData)

  use physicaldata, ONLY : facevarx, facevary, facevarz

  use Grid_data, ONLY : gr_ilo, gr_ihi, gr_jlo, gr_jhi, gr_klo, gr_khi

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(in) :: faceVar, axis, blockid
  real, intent(out), dimension(:,:,:) :: cellData

  integer :: i, j, k, ii, jj, kk
  real, dimension(:,:,:), pointer :: dataPtr

  ii = 0
  jj = 0
  kk = 0

  select case (axis)

  case(IAXIS)
    dataPtr => facevarx(faceVar,:,:,:,blockid)
    ii = 1

  case(JAXIS)
    dataPtr => facevary(faceVar,:,:,:,blockid)
    jj = 1

  case(KAXIS)
    dataPtr => facevarz(faceVar,:,:,:,blockid)
    kk = 1

  end select

  do k = 1, gr_khi - gr_klo + 1 
    do j = 1, gr_jhi - gr_jlo + 1
      do i = 1, gr_ihi - gr_ilo + 1 

        cellData(i, j, k) = 0.5 * ( dataPtr(gr_ilo+i+ii-1,gr_jlo+j+jj-1,gr_klo+k+kk-1) + dataPtr(gr_ilo+i-1,gr_jlo+j-1,gr_klo+k-1) ) 

      end do
    end do
  end do

  return

end subroutine Grid_colocateFaceData
