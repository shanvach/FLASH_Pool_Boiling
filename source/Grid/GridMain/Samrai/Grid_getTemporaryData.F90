!!****if* source/Grid/GridMain/Samrai/Grid_getTemporaryData
!!
!! NAME
!!  Grid_getTemporaryData
!!
!! SYNOPSIS
!!  #include "Flash.h"
!!
!!  Grid_getTemporaryData(integer intent(in) axis,
!!                         integer intent(in) blockId, 
!!                         real(:) intent(out) fluxes,
!!                         real(:) intent(out) area,
!!                         real(:) intent(out) dtdx,
!!                         real(:) intent(out) grav,
!!                         real(:) intent(out) ngrav,
!!                         real(:) intent(out) fict)  
!! DESCRIPTION 
!!  !!DEV: need description of this routine
!! 
!! CREATED 
!!  05/20/04, by AD
!!***

subroutine Grid_getTemporaryData(axis, blockId, size, &
                                  fluxes, area, dtdx, grav, ngrav, fict)

implicit none
  integer, intent(in) :: axis, blockId
  integer, intent(in), dimension(3) :: size
!!LBR changed dimension for compilation to 1, was :
  real, intent(out), dimension(1,size(1),size(2),size(3)) :: fluxes
  real, intent(out), dimension(size(1),size(2),size(3)) :: &
                    area,dtdx,grav,ngrav,fict

  return
end subroutine Grid_getTemporaryData





