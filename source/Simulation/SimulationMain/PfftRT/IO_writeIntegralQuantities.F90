!!****if* source/Simulation/SimulationMain/PfftRT/IO_writeIntegralQuantities
!!
!!
!!  NAME
!!    IO_writeIntegralQuantities
!!
!!  SYNOPSIS
!!    call IO_writeIntegralQuantities() 
!!                                    integer(in) :: isFirst,
!!                                    real(in)    :: simTime)
!!
!!  DESCRIPTION
!!
!!   Compute the total area of the selected countours
!!   and write them to an ASCII file.  
!!   isFirst - if 1 then write header info plus data, otherwise just write data
!!   simTime - simulation time (bogus in this application)
!!
!!
!!***

!!REORDER(4):solnData

subroutine IO_writeIntegralQuantities ( isFirst, simTime)

  use IO_data, ONLY : io_restart, io_statsFileName
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
    Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_getSingleCellVol, &
    Grid_releaseBlkPtr, Grid_fillGuardCells, Grid_getMinCellSize
  use Simulation_data
  use ut_contourSurfaceInterface, ONLY: ut_contourSurfaceAreaBlock

   use IO_data, ONLY : io_globalMe
  implicit none

#include "Flash_mpi.h"
#include "constants.h"
#include "Flash.h"
  
  
  real, intent(in) :: simTime

  integer, intent(in) :: isFirst

  integer :: lb, count, var_idx
  
  integer :: funit = 99
  integer :: error
  
  character (len=MAX_STRING_LENGTH), save :: fname 
  
  integer :: blockList(MAXBLOCKS)

  integer :: blkLimits(HIGH, MDIM), blkLimitsGC(HIGH, MDIM)

  integer, parameter ::  nGlobalSum = 3          ! Number of globally-summed quantities
  real :: gsum(nGlobalSum) !Global summed quantities
  real :: lsum(nGlobalSum) !Global summed quantities
  real :: mincellsize

  real, DIMENSION(:,:,:,:), POINTER :: solnData

  integer,parameter :: nlevels = 3
  real, dimension(nlevels) :: blkAreas
  integer :: point(MDIM)
 
  logical :: gcMask(NUNK_VARS)

  var_idx = FLSM_MSCALAR

  !**********************
  ! Fill the guardcells
  !**********************

  gcMask(:) = .false.
  gcMask(FLSM_MSCALAR) = .true.


  call Grid_fillGuardCells(CENTER, ALLDIR, &
                           maskSize=NUNK_VARS, mask=gcMask)


  if (isfirst /= 0) then
 
     var_idx = FLAM_MSCALAR

  endif

  ! Sum quantities over all locally held leaf-node blocks.
  gsum  = 0.
  lsum = 0.
  
  call Grid_getListOfBlocks(LEAF, blockList, count)
  
  do lb = 1, count
     !get the index limits of the block
     call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)

     call ut_contourSurfaceAreaBlock(nlevels,isolevels,solnData(var_idx,:,:,:), &
                                     blkLimits,blockList(lb),blkAreas)
     lsum(1) = lsum(1) + blkAreas(1)
     lsum(2) = lsum(2) + blkAreas(2)
     lsum(3) = lsum(3) + blkAreas(3)


     call Grid_releaseBlkPtr(blockList(lb), solnData)

  enddo
  

  
  ! Now the MASTER_PE sums the local contributions from all of
  ! the processors and writes the total to a file.
  
  call MPI_Reduce (lsum, gsum, nGlobalSum, FLASH_REAL, MPI_SUM, & 
       &                MASTER_PE, MPI_COMM_WORLD, error)
  
  call Grid_getMinCellSize(mincellsize)

  if (io_globalMe  == MASTER_PE) then
     
     if (isfirst == 0) then
        open (funit, file=trim(io_statsFileName), position='APPEND')
        write (funit, 12) simTime, gsum(1:3) ! Write the global sums to the file.
     else
      open (funit, file=trim(io_statsFileName))
      write (funit, 14) '# simTime', simTime 
      write (funit, 10)               &
           '#   dx                    ', &
           '    surface area flsm_1   ', &
           '    surface area flsm_2   ', &
           '    surface area flsm_3   '
      write (funit, 16) '##',mincellsize, gsum(1:3) ! Write the global sums to the file.
      write (funit, 10)               &
           '#Convolution radius       ', &
           'surface area flsm_1       ', &
           'surface area flsm_2       ', &
           'surface area flsm_3       '
     endif

14    format (2x,a9,:, 1X, es25.18)
10    format (2x,50(a25, :, 1X))
12    format (1x, 50(es25.18, :, 1x))
16    format (2x, a2,:,1x,50(es25.18, :, 1x)) 
     close (funit)          ! Close the file.
     
  endif
  
  call MPI_Barrier (MPI_Comm_World, error)
  
  !=============================================================================
  return
end subroutine IO_writeIntegralQuantities



