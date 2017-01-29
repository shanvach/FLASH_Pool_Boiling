
! Sets blockCount_start_idx from superlu_common

#include "Flash.h"

subroutine gr_slugetSubsetStartidx(blockCount_set) 

  use superlu_common, only : blockCount_start_idx, m_glob,n_glob

  use Grid_data, only :  gr_meshComm

  implicit none
!#include "constants.h"
#include "Flash_mpi.h"  
  integer, intent(in) :: blockCount_set

  ! Local variables
  integer :: aux, ierr

  ! Global numeration, by processor and block in the level:
  ! Exchange information on how many blocks have unknowns in
  ! global numeration before mine:
  call mpi_scan(blockCount_set,aux,CONSTANT_ONE,FLASH_INTEGER,MPI_SUM,gr_meshComm,ierr)

  ! blockCount_start_idx + 1 start of global numeration of MyPE subset blocks.
  blockCount_start_idx = aux - blockCount_set

  ! All reduce for block counts:
  call mpi_allreduce(blockCount_set, m_glob, CONSTANT_ONE, FLASH_INTEGER, &
                     MPI_SUM, gr_meshComm, ierr )
  m_glob = (NXB*NYB*NZB)*m_glob
  n_glob = m_glob

  return

end subroutine gr_slugetSubsetStartidx
