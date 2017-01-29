!!****if* source/Simulation/SimulationMain/unitTest/BHTree/Simulation_Perr
!!
!! NAME
!!
!!  Wind
!!
!! SYNOPSIS
!!
!!  Wind()
!!
!! DESCRIPTION
!!
!!   Calculates error of the gravitational potential.
!!
!! ARGUMENTS
!!
!! NOTES
!!   
!!
!!***

subroutine Simulation_Perr()
 
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, &
     Grid_getDeltas, Grid_releaseBlkPtr, Grid_getCellCoords, &
     Grid_getDeltas, Grid_getSingleCellVol, Grid_getListOfBlocks


  use Logfile_interface, ONLY : Logfile_stampMessage, Logfile_stamp
  use Simulation_data, ONLY : sim_MyPE, sim_Comm


 
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
#include "Eos.h"

  integer :: i, j, k, istat, ierr
  integer :: blockCount, blockID
  integer :: blockList(MAXBLOCKS)
  real :: loc_PAnlMin, loc_PNumMin, glb_PAnlMin, glb_PNumMin
  real :: dPhi, PerrMax, glb_PerrMax

  integer, dimension(2,MDIM)   :: blkLimits,blkLimitsGC
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  logical :: gcell = .true.
  character(len=MAX_STRING_LENGTH), dimension(2,2) :: strBuff


  ! Find the minimum of the analytical potential, and the value 
  ! of the numerical potential at the same point
  call Grid_getListOfBlocks(LEAF,blockList,blockCount) ! get list of LEAF blocks
  loc_PAnlMin = 0.0
  loc_PNumMin = 0.0

  do blockID = 1, blockCount
     call Grid_getBlkIndexLimits(blockList(blockID),blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blockList(blockID),solnData,CENTER)
     
     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
       do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
         do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           if (solnData(PANL_VAR, i,j,k) < loc_PAnlMin) then
             loc_PAnlMin = solnData(PANL_VAR, i,j,k)
             loc_PNumMin = solnData(GPOT_VAR, i,j,k)
           endif
         enddo
       enddo
     enddo
   enddo

   call MPI_AllReduce(loc_PAnlMin, glb_PAnlMin, 1, FLASH_REAL &
   , MPI_MIN, sim_Comm, ierr)
   call MPI_AllReduce(loc_PNumMin, glb_PNumMin, 1, FLASH_REAL &
   , MPI_MIN, sim_Comm, ierr)
   dPhi = glb_PAnlMin - glb_PNumMin

  ! Correct analytical potential (to the same minimum)
  ! and calculate the error
  PerrMax = 0.0
  do blockID = 1, blockCount
     call Grid_getBlkIndexLimits(blockList(blockID),blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blockList(blockID),solnData,CENTER)
     
     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
       do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
         do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           solnData(PANL_VAR, i,j,k) = solnData(PANL_VAR, i,j,k) - dPhi
           solnData(PERR_VAR, i,j,k) = solnData(GPOT_VAR, i,j,k) - solnData(PANL_VAR, i,j,k)
           if (abs(solnData(PERR_VAR, i,j,k)) > PerrMax) then
             PerrMax = abs(solnData(PERR_VAR, i,j,k))
           endif
         enddo
       enddo
     enddo
   enddo


   call MPI_AllReduce(PerrMax, glb_PerrMax, 1, FLASH_REAL &
   , MPI_MAX, sim_Comm, ierr)


  if ((sim_MyPE == MASTER_PE) ) then

     write (strBuff(1,1), "(A)") "analytical GPOT correction const "
     write (strBuff(1,2), "(1PE14.6)") dPhi
     call Logfile_stamp(strBuff(1:1,:), 1, 2, "[BHTree test]")

     write (strBuff(1,1), "(A)") "Max absolute GPOT error "
     write (strBuff(1,2), "(1PE14.6)") glb_PerrMax
     write (strBuff(2,1), "(A)") "Max relative GPOT error "
     write (strBuff(2,2), "(1PE14.6)") glb_PerrMax / abs(glb_PAnlMin)
     call Logfile_stamp(strBuff(1:2,:), 2, 2, "[BHTree test]")
  endif


  return
end subroutine Simulation_Perr

