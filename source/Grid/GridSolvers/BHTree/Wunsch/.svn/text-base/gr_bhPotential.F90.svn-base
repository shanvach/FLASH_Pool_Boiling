!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhPotential
!!
!! NAME
!!
!!  gr_bhPotential
!!
!!
!! SYNOPSIS
!!
!!   call gr_bhPotential(
!!          integer(in) :: idensvar,
!!          integer(in) :: ipotvar
!!          )
!!
!! DESCRIPTION
!!
!!   Computes the gravitational potential. Calls gr_bhPotentialBlock
!!   for all LEAF blocks.
!!
!! ARGUMENTS
!!
!!  idensvar - number of grid varible with density (recently not used)
!!  ipotvar  - number of grid varible with grav. potential
!!
!!***



subroutine gr_bhPotential(idensvar, ipotvar)

  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
      Grid_getListOfBlocks, Grid_updateRefinement
  use gr_bhData, ONLY : gr_bhTreeDcount, gr_bhTreeZones, &
    gr_bhTreeMyPE, gr_bhComm, gr_bhTreeBS
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Logfile_interface, ONLY : Logfile_stamp
  use gr_bhInterface, ONLY : gr_bhBuildTreeBlock
      
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer, intent(in) :: ipotvar, idensvar

  integer :: blockID, blockCount
  integer :: blockList(MAXBLOCKS)
  
  integer :: tot_zones, ierr
  real    :: tot_dcount(1:5)
  character(len=256) :: strBuff

  call Timers_start("potential")

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  gr_bhTreeDcount = 0.0
  gr_bhTreeZones = 0
  do blockID = 1, blockCount
    call gr_bhPotentialBlock(blockList(blockID), idensvar, ipotvar)
    gr_bhTreeZones = gr_bhTreeZones + gr_bhTreeBS*gr_bhTreeBS*gr_bhTreeBS
  enddo

  call MPI_Reduce(gr_bhTreeDcount,tot_dcount,5,FLASH_REAL,MPI_SUM,MASTER_PE,gr_bhComm, ierr)
  call MPI_Reduce(gr_bhTreeZones,tot_zones,1,MPI_INTEGER,MPI_SUM,MASTER_PE,gr_bhComm, ierr)

  if (gr_bhTreeMyPE == MASTER_PE) then
     write (strBuff, '("cell-cell distances: ", e10.3, ", per zone: ", f8.1)') &
     & tot_dcount(1), tot_dcount(1)/tot_zones
     call Logfile_stamp( strBuff, "[BHTree]")
     write (strBuff, '("cell-node distances: ", e10.3, ", per zone: ", f8.1)') &
     & tot_dcount(2), tot_dcount(2)/tot_zones
     call Logfile_stamp( strBuff, "[BHTree]")
     write (strBuff, '("cell-block distances: ", e10.3, ", per zone: ", f8.1)') &
     & tot_dcount(3), tot_dcount(3)/tot_zones
     call Logfile_stamp( strBuff, "[BHTree]")
     write (strBuff, '("IL cell-cell distances: ", e10.3, ", per zone: ", f8.1)') &
     & tot_dcount(4), tot_dcount(4)/tot_zones
     call Logfile_stamp( strBuff, "[BHTree]")
     write (strBuff, '("IL cell-node distances: ", e10.3, ", per zone: ", f8.1)') &
     & tot_dcount(5), tot_dcount(5)/tot_zones
     call Logfile_stamp( strBuff, "[BHTree]")
  end if

  call Timers_stop("potential")

  return
end subroutine gr_bhPotential


