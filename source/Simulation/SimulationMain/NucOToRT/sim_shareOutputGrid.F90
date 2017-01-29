#include "constants.h"
#include "Flash.h"

subroutine sim_shareOutputGrid()
  use sim_outputGridData
  use Simulation_data, ONLY : sim_meshMe
  implicit none
#include "Flash_mpi.h"


  integer :: ierr
  real,allocatable :: tmpOutGridData(:,:,:,:) !size is (0:nIOg+1, 0:nJOg+1, 0:nKOg+1, nvarsOg)
  integer,allocatable :: tmpOutGridNumCount(:,:,:) !size is (0:nIOg, 0:nJOg+1, 0:nKOg+1)
  real,allocatable :: tmpOutGridMappedVol(:,:,:) !size is (0:nIOg, 0:nJOg+1, 0:nKOg+1)

  allocate(tmpOutGridNumCount(iOgB0:iOgE0, jOgB0:jOgE0, kOgB0:kOgE0))
  call MPI_Reduce(outGridNumCount(iOgB0,jOgB0,kOgB0),tmpOutGridNumCount(iOgB0,jOgB0,kOgB0),size(outGridNumCount),FLASH_INTEGER,MPI_SUM,MASTER_PE,MPI_COMM_WORLD, ierr)
  if (sim_meshMe==MASTER_PE) then
     outGridNumCount(:,:,:) = tmpOutGridNumCount(:,:,:)
  end if
  deallocate(tmpOutGridNumCount)

  allocate(tmpOutGridMappedVol(iOgB0:iOgE0, jOgB0:jOgE0, kOgB0:kOgE0))
  call MPI_Reduce(outGridMappedVol(iOgB0,jOgB0,kOgB0),tmpOutGridMappedVol(iOgB0,jOgB0,kOgB0),size(outGridMappedVol),FLASH_REAL,MPI_SUM,MASTER_PE,MPI_COMM_WORLD, ierr)
  if (sim_meshMe==MASTER_PE) then
     outGridMappedVol(:,:,:) = tmpOutGridMappedVol(:,:,:)
  end if
  deallocate(tmpOutGridMappedVol)

  allocate(tmpOutGridData(iOgB0:iOgE0, jOgB0:jOgE0, kOgB0:kOgE0, nvarsOg))
  call MPI_Reduce(outGridData(iOgB0,jOgB0,kOgB0,1),tmpOutGridData(iOgB0,jOgB0,kOgB0,1),size(outGridData),FLASH_REAL,MPI_SUM,MASTER_PE,MPI_COMM_WORLD, ierr)
  if (sim_meshMe==MASTER_PE) then
     outGridData(:,:,:,:) = tmpOutGridData(:,:,:,:)
  end if
  deallocate(tmpOutGridData)

end subroutine sim_shareOutputGrid
