!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhExchangeTrees
!!
!! NAME
!!
!!  gr_bhExchangeTrees
!!
!!
!! SYNOPSIS
!!
!!   gr_bhExchangeTrees()
!!
!! DESCRIPTION
!!   Determines levels up to which individual block trees need to be sent to
!!   different CPUs. For a given CPU, copies all trees up to appropriate levels
!!   to a single message and sends it to the CPU.
!!   
!!
!! ARGUMENTS
!!
!!
!!***



subroutine gr_bhExchangeTrees()

  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
      Grid_getListOfBlocks, Grid_updateRefinement, &
      Grid_getBlkBoundBox
  use gr_bhData, ONLY : GR_TREE_BND_PERIODIC, &
    gr_bhLx, gr_bhLy, gr_bhLz, p_message, gr_bhTreeLevels, &
    gr_bhLocSentTreeLevels, gr_bhTreeNodetype, gr_bhTreeBCen, &
    gr_bhTreeArray, gr_bhComm, gr_bhTreeBS, gr_bhTreeLimAngle2, &
    gr_bhTreeMyPE, gr_bhTreeNumProcs, gr_bhBndType, &
    gr_bhTreeLrefine, gr_bhTreeCellSize, gr_bhTreeLnblocks, &
    gr_bhTreeDiag2, gr_bhLocRecvTreeLevels

  use Timers_interface, ONLY : Timers_start, Timers_stop
  use gr_bhInterface, ONLY : gr_bhGetTreeSize
  use gr_bhLocalInterface, ONLY : gr_bhGetTreePos
  use tree, ONLY : nodetype, lrefine, child, mchild, lnblocks
  use Logfile_interface, ONLY : Logfile_stamp

      
  implicit none
#include "constants.h"
#include "Flash_mpi.h"

  integer :: blockID, blockCount
  integer :: blockList(MAXBLOCKS)
  integer :: ierr, istat, i, j, k, ii
  integer :: lBID, rBID, l, ts
  real    :: ldiag, rdiag, dist2, mindist2, dx, dy, dz
  integer :: ms_size(0:gr_bhTreeNumProcs), mr_size(0:gr_bhTreeNumProcs)
  integer :: stats(MPI_STATUS_SIZE,2*(gr_bhTreeNumProcs-1)), reqs(2*(gr_bhTreeNumProcs-1))
  integer :: lcor, rcor, ls1, ls2, ls3, rs1, rs2, rs3
  integer :: multi(1:gr_bhTreeLevels), level, pos
  type(p_message), save, allocatable :: messages_send(:), messages_recv(:)

  call Timers_start("exchange_tree")


  ! determine levels of trees to be sent to different CPUs
  ! find distance of block centers minus both block space diagonals
  call Timers_start("et: det_gr_bhTreeLevels")
  do i = 0,gr_bhTreeNumProcs-1
    ms_size(i) = 0
    if (i /= gr_bhTreeMyPE) then
      ! find the minimum distance of the local block to any block on CPU i
      do lBID = 1, gr_bhTreeLnblocks(gr_bhTreeMyPE) ! local blocks
        if (nodetype(lBID) .ne. 1) cycle
        mindist2 = 1d99
        do rBID = 1, gr_bhTreeLnblocks(i) ! remote blocks (on CPU i)
          if (gr_bhTreeNodetype(rBID, i) .ne. 1) cycle

          dx = abs(gr_bhTreeBCen(IAXIS,rBID,i)-gr_bhTreeBCen(IAXIS,lBID,gr_bhTreeMyPE)) 
          if (gr_bhBndType(1) .EQ. GR_TREE_BND_PERIODIC) dx = min(abs(dx), abs(dx+gr_bhLx), abs(dx-gr_bhLx)) 
          dx = dx &
             - 0.5*gr_bhTreeBS*gr_bhTreeCellSize(gr_bhTreeLrefine(rBID, i), IAXIS) &
             - 0.5*gr_bhTreeBS*gr_bhTreeCellSize(gr_bhTreeLrefine(lBID, gr_bhTreeMyPE), IAXIS)
          dy = abs(gr_bhTreeBCen(JAXIS,rBID,i)-gr_bhTreeBCen(JAXIS,lBID,gr_bhTreeMyPE)) 
          if (gr_bhBndType(3) .EQ. GR_TREE_BND_PERIODIC) dy = min(abs(dy), abs(dy+gr_bhLy), abs(dy-gr_bhLy))
          dy = dy &
             - 0.5*gr_bhTreeBS*gr_bhTreeCellSize(gr_bhTreeLrefine(rBID, i), JAXIS) &
             - 0.5*gr_bhTreeBS*gr_bhTreeCellSize(gr_bhTreeLrefine(lBID, gr_bhTreeMyPE), JAXIS)
          dz = abs(gr_bhTreeBCen(KAXIS,rBID,i)-gr_bhTreeBCen(KAXIS,lBID,gr_bhTreeMyPE)) 
          if (gr_bhBndType(5) .EQ. GR_TREE_BND_PERIODIC) dz = min(abs(dz), abs(dz+gr_bhLz), abs(dz-gr_bhLz))
          dz = dz &
             - 0.5*gr_bhTreeBS*gr_bhTreeCellSize(gr_bhTreeLrefine(rBID, i), KAXIS) &
             - 0.5*gr_bhTreeBS*gr_bhTreeCellSize(gr_bhTreeLrefine(lBID, gr_bhTreeMyPE), KAXIS)
          dx = max(dx, 0.0)
          dy = max(dy, 0.0)
          dz = max(dz, 0.0)

          dist2 = dx*dx + dy*dy + dz*dz
          if (dist2 .lt. mindist2) mindist2 = dist2
        enddo

        
        ! set the level of local tree lBID sent to CPU i
        gr_bhLocSentTreeLevels(lBID, i) = gr_bhTreeLevels
        if (mindist2 > 1.d-99) then
          do l = 0,gr_bhTreeLevels
            if (gr_bhTreeDiag2(l+gr_bhTreeLrefine(lBID,gr_bhTreeMyPE))/(mindist2+1.d-99) < gr_bhTreeLimAngle2) then
              gr_bhLocSentTreeLevels(lBID, i) = l
              exit
            endif
          enddo
        endif

        ! add tree size to the ms_size
        ms_size(i) = ms_size(i) + gr_bhGetTreeSize(gr_bhLocSentTreeLevels(lBID, i))
      enddo
    endif
  enddo
  call Timers_stop("et: det_gr_bhTreeLevels")

  call Timers_start("et: comm_gr_bhTreeLevels")

  if (gr_bhTreeMyPE == MASTER_PE) then
     call Logfile_stamp( "ET before STL", "[BHTree]")
  end if
  ! sent sentTreeLevels to recvTreeLevels :)
  j = 1 ! tag of the message, use also for indexing reqs and stats arrays
  do i = 0,gr_bhTreeNumProcs-1
    if (i /= gr_bhTreeMyPE) then
      call MPI_Irecv(gr_bhLocRecvTreeLevels(1,i), MAXBLOCKS, MPI_INTEGER, i, 1, gr_bhComm, reqs(j), ierr)
      call MPI_Isend(gr_bhLocSentTreeLevels(1,i), MAXBLOCKS, MPI_INTEGER, i, 1, gr_bhComm, reqs(gr_bhTreeNumProcs-1+j), ierr)
      j = j + 1
    endif
  enddo
  call MPI_WaitAll(2*(gr_bhTreeNumProcs-1), reqs, stats, ierr)
  

  if (gr_bhTreeMyPE == MASTER_PE) then
     call Logfile_stamp( "ET after STL", "[BHTree]")
  end if

  call Timers_stop("et: comm_gr_bhTreeLevels")
  call Timers_start("et: comm_trees")

  ! build messages to be sent, fill them with local trees
  allocate(messages_send(0:gr_bhTreeNumProcs), stat=istat)
  do i = 0,gr_bhTreeNumProcs-1
    if (i /= gr_bhTreeMyPE) then
      ! create message to be sent to i-th CPU
      nullify(messages_send(i)%p)
      allocate(messages_send(i)%p(1:ms_size(i)), stat=istat)
      j = 1
      do lBID = 1, gr_bhTreeLnblocks(gr_bhTreeMyPE)
        if (nodetype(lBID) .ne. 1) cycle
        do k = 1, gr_bhGetTreeSize(gr_bhLocSentTreeLevels(lBID, i))
          messages_send(i)%p(j) = gr_bhTreeArray(gr_bhTreeMyPE, lBID)%p(k)
          j = j + 1
        enddo
      enddo
    endif
  enddo

  
  ! prepare space for messages to be recieved
  allocate(messages_recv(0:gr_bhTreeNumProcs), stat=istat)
  do i = 0,gr_bhTreeNumProcs-1
    if (i /= gr_bhTreeMyPE) then
      j = 0
      do rBID = 1, gr_bhTreeLnblocks(i)
        if (gr_bhTreeNodetype(rBID, i) .ne. 1) cycle
        j = j + gr_bhGetTreeSize(gr_bhLocRecvTreeLevels(rBID, i))
      enddo
      mr_size(i) = j
      nullify(messages_recv(i)%p)
      allocate(messages_recv(i)%p(1:mr_size(i)))
    endif
  enddo

  ! to ensure buffers for messages are allocated
  ! may not be necessary
  call MPI_Barrier(gr_bhComm, ierr)
  if (gr_bhTreeMyPE == MASTER_PE) then
     call Logfile_stamp( "ET tree mes allocated", "[BHTree]")
  end if
  
  ! exchange messages
  j = 1 ! tag of the message, use also for indexing reqs and stats arrays
  do i = 0,gr_bhTreeNumProcs-1
    if (i /= gr_bhTreeMyPE) then
      call MPI_Irecv(messages_recv(i)%p, mr_size(i), FLASH_REAL, i, 1, gr_bhComm, reqs(j), ierr)
      call MPI_Isend(messages_send(i)%p, ms_size(i), FLASH_REAL, i, 1, gr_bhComm, reqs(gr_bhTreeNumProcs-1+j), ierr)
      j = j + 1
    endif
  enddo
  call MPI_WaitAll(2*(gr_bhTreeNumProcs-1), reqs, stats, ierr)
  
  if (gr_bhTreeMyPE == MASTER_PE) then
     call Logfile_stamp( "ET tree exchanged", "[BHTree]")
  end if
  
  ! copy data from messages to trees
  do i = 0,gr_bhTreeNumProcs-1
    if (i /= gr_bhTreeMyPE) then
      j = 1
      do rBID = 1, gr_bhTreeLnblocks(i)
        if (gr_bhTreeNodetype(rBID, i) .ne. 1) cycle
        ts = gr_bhGetTreeSize(gr_bhLocRecvTreeLevels(rBID, i))
        nullify(gr_bhTreeArray(i,rBID)%p)
        allocate(gr_bhTreeArray(i,rBID)%p(1:ts), stat=istat)
        if (istat /= 0) call Driver_abortFlash("could not allocate tree in gr_bhExchangeTrees.F90")
        do k = 1, ts
          gr_bhTreeArray(i, rBID)%p(k) = messages_recv(i)%p(j)
          j = j + 1
        enddo
      enddo
    endif
  enddo

  ! destroy messages
  do i = 0,gr_bhTreeNumProcs-1
    if (i /= gr_bhTreeMyPE) then
      deallocate(messages_send(i)%p)
      deallocate(messages_recv(i)%p)
    endif
  enddo
  deallocate(messages_send)
  deallocate(messages_recv)
  call Timers_stop("et: comm_trees")
  
  call Timers_stop("exchange_tree")

  return
end subroutine gr_bhExchangeTrees


