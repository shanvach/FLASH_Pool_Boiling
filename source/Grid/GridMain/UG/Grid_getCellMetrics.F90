!!****if* source/Grid/GridMain/UG/UGStretched/Grid_getCellMetrics
!!
!! NAME
!!  Grid_getCellMetrics
!!
!! SYNOPSIS
!!
!!  Grid_getCellMetrics(integer(IN) :: axis,
!!		        integer(IN) :: blockId,
!!                      integer(IN) :: edge,
!!                      logical(IN) :: guardcell,
!!                      real(OUT)   :: metrics(size),
!!                      integer(IN) :: size)
!!
!! DESCRIPTION
!!
!!    This subroutine is an accessor function that gets the derivative 
!!    metric coefficients of the cells in a given block.  Metrics are 
!!    retrieved one axis and location within a cell at a time, 
!!    meaning you can get the i, j, or k derivative direction metrics
!!    at the eastern cell face with one call. If you want all the metrics, 
!!    for all the axes, and all the cell locations, you need to call 
!!    Grid_getCellMetrics 12 times, one for each axis and location.
!!    The code carries metrics for each axis at cell centers as well as faces.
!!    It is possible to get metrics for CENTER, FACEX, FACEY, or FACEZ.
!!
!!
!!
!!
!! ARGUMENTS
!!            
!!   axis - specifies the integer index direction of the metrics being retrieved.
!!          axis can have one of three different values, IAXIS, JAXIS or KAXIS 
!!          (defined in constants.h as 1,2 and 3)
!!
!!   blockId - integer block number
!!
!!   edge - integer value specifing the location within the cells to retrieve the 
!!          associated metrics,
!!          CENTER cell centered metrics
!!          FACEX  face centered metrics on faces along IAXIS
!!          FACEY  face centered metrics on faces along JAXIS
!!          FACEZ  face centered metrics on faces along IAXIS 
!!
!!   guardcell - logical input. If true metrics for guardcells are returned
!!          along with the interior cells, if false, only the interior metrics 
!!          are returned.
!!
!!          
!!   metrics - The array holding the data returning the metric values
!!          metrics must be at least as big as "size" (see below)
!!           
!!   size - integer specifying the size of the metrics array.
!!          If guardcell true then size =  interior cells + 2*guardcells
!!          otherwise size = number of interior cells
!!
!!
!!  EXAMPLE 
!!
!!  1. Getting cell centered metrics in the x direction (dx @ cell centers)
!!
!!   #include "constants.h"
!!   #include "Flash.h"
!!      
!!      integer :: size
!!      integer :: metrics(size) !sized to be number of metrics returned
!!      
!!      do i=1, localNumBlocks
!!
!!          !get the index limits of the block
!!          call Grid_getBlkIndexLimits(i, blkLimits, blkLimitsGC)
!!
!!          !holds the number of cells returned in idir
!!          coordSize = blkLimitsGC(HIGH, IAXIS)
!!          call Grid_getCellMetrics(IAXIS, i, CENTER, .true., metrics, size) 
!!
!!     end do    
!!
!!  2. Getting northern face centered metrics in the y direction (dy @ cell north faces)
!! 
!!   #include "constants.h"
!!   #include "Flash.h"
!!
!!      
!!      integer :: size
!!      integer :: metrics(size) !sized to be number of metrics returned
!!      
!!      do i=1, localNumBlocks
!!
!!          !get the index limits of the block
!!          call Grid_getBlkIndexLimits(i, blkLimits, blkLimitsGC)
!!
!!          !holds the number of cells returned in ydir
!!          size = blkLimitsGC(HIGH, JAXIS)
!!          call Grid_getCellMetrics(JAXIS, i, FACEZ, .true., metrics, size) 
!!
!!     end do    
!!
!!
!!***

subroutine Grid_getCellMetrics(axis, blockId, edge, guardcell, metrics, size)

  implicit none
  integer, intent(in) :: axis
  integer, intent(in) :: blockId, edge
  integer, intent(in) :: size
  logical, intent(in) :: guardcell
  real,    intent(out), dimension(size) :: metrics

  print*,"Get Cell Metrics : invalid operation for uniform grid "
  call Driver_abortFlash("Get Cell Metrics : invalid operation for uniform grid ")  

  return
end subroutine Grid_getCellMetrics

