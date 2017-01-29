subroutine Cool_unitTest( fileUnit, perfect)

  use Grid_interface, ONLY: Grid_getListOfBlocks, Grid_getBlkPtr, &
       Grid_releaseBlkPtr, Grid_getLocalNumBlks, Grid_getBlkIndexLimits
use Cool_data, ONLY : cl_meshMe
  use Cool_interface, ONLY : Cool

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer, intent(IN) ::  fileUnit
  logical, intent(INOUT) :: perfect

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: blockCount, lb, blockID, localNumBlocks
  integer :: blockList(MAXBLOCKS), i, j, k, ierr

  real,pointer,dimension(:,:,:,:) :: solnData

  real :: dt, time, dens, temp, dedt

  call Grid_getLocalNumBlks(localNumBlocks)
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  ! Set the timestep to something fiducial for cosmological/cluster
  ! situations.

  dt = 3.1557e14 ! 10 Myr
  time = 0.0

  call Cool(blockCount, blockList, dt, time)

  if (cl_meshMe == 0) print *, "Done calling the cooling function"

  ! Now dump crap to a text file

  if (cl_meshMe == 0) then

     open (9, file="cooling_fctn.dat")

     do lb = 1, blockCount

        blockID = blockList(lb)

        call Grid_getBlkPtr(blockID,solnData)
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                 
                 dens = solnData(DENS_VAR,i,j,k)
                 temp = solnData(TEMP_VAR,i,j,k)
                 dedt = solnData(COOL_VAR,i,j,k)
                 
                 write(9,*) dens, temp, dedt

              enddo
           enddo
        enddo
        
        call Grid_releaseBlkPtr(blockID, solnData)

     enddo

  endif

  call mpi_barrier(MPI_COMM_WORLD, ierr)

  return

end subroutine Cool_unitTest
