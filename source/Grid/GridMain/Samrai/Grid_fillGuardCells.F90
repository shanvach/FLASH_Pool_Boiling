!!****if* source/Grid/GridMain/Samrai/Grid_fillGuardCells
!!
!! NAME
!!  Grid_fillGuardCells
!!
!! SYNOPSIS
!!  #include Flash.h
!!
!!  Grid_fillGuardCells(integer(IN) :: myPE, 
!!                          integer(IN) :: gridDataStruct,
!!                          integer(IN) :: idir,
!!                 optional,integer(IN) :: minLayers,
!!                  logical,optional(IN):: mask,
!!                 optional,integer(IN) :: level)
!!  
!! DESCRIPTION 
!!  
!!  Fills guard cells for all finest blocks
!!  DOC: needs more or a description
!!  
!! ARGUMENTS 
!!  
!!  myPE         - local processor number
!!  gridDataStruct - some Grid implementation have multiple data structure
!!                   this argument indicates which one to fill 
!!  idir - direction of guardcell fill.  User can specify ALLDIR for all (x,y,z)
!!         directions, or if for example the algorithm only does one directional
!!         sweep at a time then time can be saved by filling only the guardcell
!!         direction that is needed
!!  minLayers - number of guardcell layers requested for all directions.
!!              DEV: minlayers currently ignored in the Samrai implementation.
!!  mask - If this argument is present, then selected variables have their
!!         guard cells filled. mask has same indices as the grid data strucures
!!         and if a variable's guard cells are to be filled, its corresponding 
!!         element in mask is true
!!  level - This optional argument is for future use.
!!
!!
!! NOTES
!!
!!  DEV: at some point we may want to add in SimTime as a parameter
!!  for problems with time dependent boundary conditions
!!  
!!  DEV: verify interface for filling guardcells at only one level.
!!  not yet implemented, may want additional or optional arguments like
!!  validDataFinestOnly - integer if true, fill GC using data from
!!                        finest level neighbors only.
!!                        if false, fill GC using data
!!                        from neighbors on the same level,
!!                        even if they are covered by finer blocks.
!!***

subroutine Grid_fillGuardCells(myPE, gridDataStruct, idir,minLayers, mask, level)
  use Grid_data, ONLY: gr_convertToConsvdForMeshCalls
  implicit none
#include "constants.h"
#include "Flash.h"
  integer, intent(in) :: myPE
  integer, intent(in) :: gridDataStruct
  integer :: mapperType = NO_MAPPER
  integer, intent(in) :: idir
  integer, optional,intent(in) :: minLayers
  logical,dimension(:),optional,intent(in) :: mask
  integer, optional,intent(in) :: level !! Has no meaning for SAMRAI
 
  if(gr_convertToConsvdForMeshCalls) mapperType = CONS_TO_PRIM
  call samrai_fill_guard_cells(myPE, gridDataStruc, mapperType, idir)

end subroutine Grid_fillGuardCells
