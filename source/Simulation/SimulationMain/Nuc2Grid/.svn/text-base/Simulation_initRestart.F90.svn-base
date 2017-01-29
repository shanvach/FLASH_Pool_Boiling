!!****if* source/Simulation/SimulationMain/Nuc2Grid/Simulation_initRestart
!!
!! NAME
!!  Simulation_initRestart
!!
!! SYNOPSIS
!!  call Simulation_initRestart(
!!
!! DESCRIPTION
!!  This is where the user should place code for a setup that needs to adjust
!!  data on a restart, particularly if grid data, grid metadata or particle 
!!  data needs to be changed on restarting.
!!
!! ARGUMENTS
!!
!!***

!!REORDER(4):solnData

#include "Flash.h"
#include "constants.h"

subroutine Simulation_initRestart()

  use Simulation_data, ONLY : sim_unkCellWeight
  use Grid_interface,      ONLY : Grid_getListOfBlocks, &
                                  Grid_getBlkIndexLimits,&
                                  Grid_getBlkPtr,&
                                  Grid_releaseBlkPtr
                                  
  implicit none

  real, DIMENSION(:,:,:,:), POINTER :: solnData

  integer :: lb, blockCount
  integer :: blockList(MAXBLOCKS)
  integer :: blkLimits(HIGH, MDIM), blkLimitsGC(HIGH, MDIM)
  integer :: ivar,i,j,k

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  do lb = 1, blockCount
     !get the index limits of the block
     call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)

     do ivar = 1,NUNK_VARS

        if(sim_unkCellWeight(ivar) .NE. 1.0) then
           do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
              do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                 do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
!!$                    if (sim_unkCellWeight(ivar) > 0.0 .AND. sim_unkCellWeight(ivar)<1.0e-8) then
!!$                       select case (ivar)
!!$                       case (C_SPEC,O_SPEC)
!!$                          solnData(ivar,i,j,k) = 0.5
!!$                       end select
!!$                    end if
                    solnData(ivar,i,j,k) = solnData(ivar,i,j,k) * sim_unkCellWeight(ivar)
                 enddo
              enddo
           enddo
        end if
     enddo
     call Grid_releaseBlkPtr(blockList(lb), solnData)

  enddo

end subroutine Simulation_initRestart
