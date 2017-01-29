!!****if* source/Grid/GridMain/Samrai/Grid_dump
!!
!! NAME
!!  Grid_dump
!!
!! SYNOPSIS
!!  #include "Flash.h"
!!
!!  call Grid_dump(integer(IN) :: var(num),
!!                 integer(IN) :: num,
!!                 integer(IN) :: blockID
!!                 logical(IN) :: gcell)
!!  
!! DESCRIPTION 
!!  
!! Binary dumps the variables specified in "var". Can be done from 
!! anywhere in the code, and is useful for diagnostic purposes
!!  
!! ARGUMENTS 
!!
!!  var :: array containing the indices of the variables to be dumped
!!  num :: number of variables being dumped.
!!  blk :: number of block to be dumped
!!  gcell :: indicates whether to include guardcells in the dump.
!!             
!! EXAMPLES 
!!  if num = 3, and var(1)=DENS_VAR, var(2)=PRES_VAR and var(3)=TEMP_VAR
!!  then the current values of density, pressure and temperature will be 
!!  dumped.
!!
!!***

subroutine Grid_dump(var, num, blockID, gcell)

  use Grid_data, ONLY : gr_maxPatches  

  implicit none

  integer, intent(IN) :: blockID, num
  integer, intent(IN) :: var
  logical, intent(IN) :: gcell


  integer :: levelNum, patchNum
  integer :: varIDSamrai
  
  !translate blockID to patchNum and levelNum
  levelNum = blockID / gr_maxPatches

  patchNum = mod(blockID, gr_maxPatches)

  varIDSamrai = varID - 1

  call samrai_dump(patchNum, levelNum, varIDSamrai, gcell)


  return
end subroutine Grid_dump
