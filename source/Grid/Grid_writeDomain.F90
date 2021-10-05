!!****if* source/Grid/GridMain/Grid_writeDomain
!!
!! NAME
!!
!!  Grid_writeDomain
!!
!!
!! SYNOPSIS
!!
!!  call Grid_writeDomain()
!!
!!
!! DESCRIPTION
!!
!! This function writes the grid information to an hdf5 file to store the 
!! Paramesh or Uniform Grid / Regular Grid cell coordinates (Left, Center, Right) 
!! and the cell metrics for later use in post-processing FLASH simulations.
!!
!! Currently only supports hdf5 IO
!!
!! This fuction is intended to be called after a IO_output function
!!
!! ARGUMENTS
!!
!!  none
!!
!!***
 
subroutine Grid_writeDomain(fileNumber)

  implicit none
  integer, optional, intent(IN) :: fileNumber

end subroutine Grid_WriteDomain
