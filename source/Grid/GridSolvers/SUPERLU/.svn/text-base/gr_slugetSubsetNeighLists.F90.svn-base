
! Sets arrays neighProc_blkcnt,neighProc_blkList,neighProc_blkList_idx
! defined in superlu_common


#include "constants.h"
#include "Flash.h"

subroutine gr_slugetSubsetNeighLists(blockList_set,blockCount_set)

  use superlu_common, only : neighProcsCount,neighProcsList,     &
                             neighProc_blkcnt,neighProc_blkList, &
                             blockCount_start_idx,neighProc_blkList_idx

  use Grid_data, only : gr_meshMe,gr_meshNumProcs,gr_meshComm

  implicit none
#include "Flash_mpi.h"
  integer, intent(in) :: blockList_set(MAXBLOCKS),blockCount_set

  ! Local Variables
  integer :: recvs(gr_meshNumProcs),sts(MPI_STATUS_SIZE,gr_meshNumProcs),vecBuff(2),TAG,ierr
  integer :: p, b, neighProc_idx, neighProc_blocks, blockID

  ! Allocate data arrays:
if (allocated(neighProc_blkcnt)) deallocate(neighProc_blkcnt)
if (allocated(neighProc_blkList)) deallocate(neighProc_blkList)
  allocate(neighProc_blkcnt(CONSTANT_TWO,neighProcsCount)) 
  allocate(neighProc_blkList(MAXBLOCKS,neighProcsCount)); neighProc_blkList = 0
  
  ! Send and receive blockCount_set among list of neighProcs => fill neighProc_blkcnt
  recvs  = MPI_REQUEST_NULL
  TAG    = CONSTANT_ONE
  if (neighProcsCount > 0) then 
     do p=1,neighProcsCount
        call MPI_Irecv(neighProc_blkcnt(1:2,p),CONSTANT_TWO,FLASH_INTEGER, &
                       neighProcsList(p),TAG,gr_meshComm,recvs(p),ierr)
     enddo
     vecBuff(CONSTANT_ONE) = blockCount_start_idx
     vecBuff(CONSTANT_TWO) = blockCount_set
     do p=1,neighProcsCount
        call MPI_Ssend(vecBuff(1:2),CONSTANT_TWO,FLASH_INTEGER,neighProcsList(p), &
                       TAG,gr_meshComm,ierr) !sends(p)
     enddo
  end if

  call MPI_WaitAll(neighProcsCount, recvs, sts, ierr)
  if (ierr /= MPI_SUCCESS) then
     call Driver_abortFlash("Send MPI_Waitall error for neighProc_blkcnt")
  endif

  ! Send and receive blockList_set among list of neighProcs  => fill neighProc_blkidx
  ! With these we find in lo
  recvs  = MPI_REQUEST_NULL
  TAG    = CONSTANT_TWO
  if (neighProcsCount > 0) then 
     do p=1,neighProcsCount
        neighProc_blocks = neighProc_blkcnt(CONSTANT_TWO,p)
        call MPI_Irecv(neighProc_blkList(1:neighProc_blocks,p),neighProc_blocks,FLASH_INTEGER,&
                       neighProcsList(p),TAG,gr_meshComm,recvs(p),ierr)
     enddo
     do p=1,neighProcsCount
        call MPI_Ssend(blockList_set(1:blockCount_set),blockCount_set,&
                       FLASH_INTEGER,neighProcsList(p),TAG,gr_meshComm,ierr) !sends(p)
     enddo
  end if

  call MPI_WaitAll(neighProcsCount, recvs, sts, ierr)
  if (ierr /= MPI_SUCCESS) then
     call Driver_abortFlash("Send MPI_Waitall error for neighProc_blkidx")
  endif

  ! Now build Local to global block indexes of gr_meshMe and neighProcs.
  ! Make it simple, stupid!!
  ! gr_meshMe is number 1 on the list
if (allocated(neighProc_blkList_idx)) deallocate(neighProc_blkList_idx)
  allocate(neighProc_blkList_idx(MAXBLOCKS,neighProcsCount)); neighProc_blkList_idx = 0
  do p=1,neighProcsCount
     neighProc_idx    = neighProc_blkcnt(1,p) 
     neighProc_blocks = neighProc_blkcnt(2,p)
     do b = 1,neighProc_blocks
        blockID = neighProc_blkList(b,p)
        neighProc_blkList_idx(blockID,p) = neighProc_idx + b
     enddo
  enddo

  return

end subroutine gr_slugetSubsetNeighLists
