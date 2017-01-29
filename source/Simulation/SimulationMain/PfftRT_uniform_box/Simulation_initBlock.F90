subroutine Simulation_initBlock(blockID)


!==============================================================================
! Set up a Gaussian in the middle of the domain, as a test problem
!==============================================================================

  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_putPointData, &
                             Grid_getCellCoords, Grid_getDomainBoundBox
  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(in) :: blockID
  
  integer :: i, j, k, d
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, allocatable, dimension(:) :: iCoords, jCoords, kCoords
  integer, dimension(MDIM) :: cell
  integer :: isize, jsize, ksize
  real, dimension(2, MDIM) :: bbox
  real, dimension(NDIM) :: domain_ctr, gwidth
  real :: dens, flam, flsm, auxv, gaus, arg1, arg2, arg3
  interface
    real function boxfn(coord, ctr, wid)
        implicit none
        real, intent(IN) :: coord, ctr, wid
    end function boxfn
  end interface

!==============================================================================
! Initialize a uniform box
!==============================================================================

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  isize = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
  jsize = blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
  ksize = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1
  allocate(iCoords(isize))
  allocate(jCoords(jsize))
  allocate(kCoords(ksize))
  call Grid_getCellCoords(IAXIS,blockID,CENTER,.false.,iCoords,isize)
  call Grid_getCellCoords(JAXIS,blockID,CENTER,.false.,jCoords,jsize)
  call Grid_getCellCoords(KAXIS,blockID,CENTER,.false.,kCoords,ksize)

  call Grid_getDomainBoundBox(bbox)
  do d = 1, NDIM

    domain_ctr(d) = 0.5 * ( bbox(2, d) + bbox(1,d) )
    gwidth(d) = gwidth_rel(d) * 0.5 * ( bbox(2, d) - bbox(1,d) )

  end do
  
  arg2 = 0.0
  arg3 = 0.0
  do i = 1, isize
    cell(1) = i
    arg1 = boxfn(iCoords(i), domain_ctr(1), gwidth(1))

    do j = 1, jsize
      cell(2) = j
      arg2 = boxfn(jCoords(j), domain_ctr(2), gwidth(2))
      
      do k = 1, ksize
        cell(3) = k
        arg3 = boxfn(kCoords(k), domain_ctr(3), gwidth(3))

        flam = arg1 * arg2 * arg3
        dens = 0.0
        flsm = 0.0
        auxv = 0.0
        gaus = 0.0
        
        call Grid_putPointData(blockId, CENTER, FLAM_VAR, INTERIOR, cell, flam)
        call Grid_putPointData(blockId, CENTER, FLSM_VAR, INTERIOR, cell, flsm)
        call Grid_putPointData(blockId, CENTER, DENS_VAR, INTERIOR, cell, dens)
        call Grid_putPointData(blockId, CENTER, AUXV_VAR, INTERIOR, cell, auxv)
        call Grid_putPointData(blockId, CENTER, GAUS_VAR, INTERIOR, cell, gaus)
        
      end do
    end do
  end do
  
  deallocate(iCoords)
  deallocate(jCoords)
  deallocate(kCoords)

end subroutine Simulation_initBlock

real function boxfn(coord, ctr, wid)

    implicit none
    real, intent(IN) :: coord, ctr, wid

    if ( abs( coord - ctr ) .le. wid ) then
        boxfn = 0.5 / wid
    else
        boxfn = 0.0
    endif

end function boxfn



