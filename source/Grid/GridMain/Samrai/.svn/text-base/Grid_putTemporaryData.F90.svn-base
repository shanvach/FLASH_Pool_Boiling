!!****if* source/Grid/GridMain/Samrai/Grid_putTemporaryData
!!
!! NAME
!!  Grid_putTemporaryData
!!
!! SYNOPSIS
!!  #include "Flash.h"
!!
!!  Grid_putTemporaryData(integer intent(in) axis,
!!                         integer intent(in) blockId, 
!!                         real(NFLUXES, size(1), size(2), size(3)) intent(in) fluxes,
!!                         real(size(1),size(2), size(3)) intent(in) area,
!!                         real(size(1),size(2), size(3)) intent(in) dtdx,
!!                         real(size(1),size(2), size(3)) intent(in) grav,
!!                         real(size(1),size(2), size(3)) intent(in) ngrav,
!!                         real(size(1),size(2), size(3)) intent(in) fict)
!!  
!! DESCRIPTION 
!! 
!! CREATED 
!!
!!***

subroutine Grid_putTemporaryData(axis, blockId, size, &
                                  fluxes, area, dtdx, grav, ngrav, fict)

implicit none
  integer, intent(in) :: axis, blockId
  integer, intent(in), dimension(3) :: size
!!LBR changed this first dimension to 1 for stub compilation -- was :
  real, intent(in), dimension(1,size(1),size(2),size(3)) :: fluxes
  real, intent(in), dimension(size(1),size(2),size(3)) :: &
                    area,dtdx,grav,ngrav,fict
  return
end subroutine Grid_putTemporaryData





