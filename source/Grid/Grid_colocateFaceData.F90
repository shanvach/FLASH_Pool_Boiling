!!****if* source/Grid/GridMain/UG/Grid_writeDomain
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
     !!      and only for the uniform grid, and derived, grids
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

  implicit none

end subroutine Grid_colocateFaceData
