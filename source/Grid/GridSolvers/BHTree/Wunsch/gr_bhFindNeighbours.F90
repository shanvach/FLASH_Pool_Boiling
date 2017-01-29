!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhFindNeighbours
!!
!! NAME
!!
!!  gr_bhFindNeighbours
!!
!!
!! SYNOPSIS
!!
!!   gr_bhFindNeighbours()
!!
!! DESCRIPTION
!!
!!   Finds neighbours (including the diagonal ones) for all blocks in the simulation.
!!   The gr_bhTreeSurbox array includes 2x27 elements for each block. They give block number
!!   and CPU of neighbour blocks. Neighbours are numbered by the following scheme:
!!
!!          k-1           k               k+1
!!    ^  07 08 09      16 17 18        25 26 27
!!    |  04 05 06      13 14 15        22 23 24
!!    j  01 02 03      10 11 12        19 20 21
!!       i ->
!!
!!   (element 14 refers to itself).
!!
!! ARGUMENTS
!!
!!
!!***




subroutine gr_bhFindNeighbours()

  use gr_bhData, ONLY: tr_surbox => gr_bhTreeSurbox, tr_neigh => gr_bhTreeNeigh, &
       gr_bhTreeNumProcs, gr_bhTreeLrefine, gr_bhTreeNodeType
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer :: block_no, cpu, i

  do cpu = 0, gr_bhTreeNumProcs-1
    do block_no = 1,MAXBLOCKS

      !print *, "STARTING WITH BLOCK: ", block_no, cpu, gr_bhTreeNodetype(block_no,cpu)
      if (gr_bhTreeNodetype(block_no,cpu) .ne. 1) cycle  ! 1 = LEAF block

      do i = 1,27
        tr_surbox(1, i, block_no, cpu) = -1
        tr_surbox(2, i, block_no, cpu) = -1
      enddo

      ! central point
      tr_surbox(1,14,block_no,cpu) = block_no
      tr_surbox(2,14,block_no,cpu) = cpu
      !print *, "CENTRAL POINT: ", block_no, cpu

      ! face neighbours: 05,11,13,15,17,23
      tr_surbox(1, 5,block_no,cpu) = tr_neigh(1,5,block_no,cpu) ! 5 = low z neigh
      tr_surbox(2, 5,block_no,cpu) = tr_neigh(2,5,block_no,cpu)
      tr_surbox(1,11,block_no,cpu) = tr_neigh(1,3,block_no,cpu) ! 3 = low y neigh
      tr_surbox(2,11,block_no,cpu) = tr_neigh(2,3,block_no,cpu)
      tr_surbox(1,13,block_no,cpu) = tr_neigh(1,1,block_no,cpu) ! 1 = low x neigh
      tr_surbox(2,13,block_no,cpu) = tr_neigh(2,1,block_no,cpu)
      tr_surbox(1,15,block_no,cpu) = tr_neigh(1,2,block_no,cpu) ! 2 = high x neigh
      tr_surbox(2,15,block_no,cpu) = tr_neigh(2,2,block_no,cpu)
      tr_surbox(1,17,block_no,cpu) = tr_neigh(1,4,block_no,cpu) ! 5 = high y neigh
      tr_surbox(2,17,block_no,cpu) = tr_neigh(2,4,block_no,cpu)
      tr_surbox(1,23,block_no,cpu) = tr_neigh(1,6,block_no,cpu) ! 5 = high z neigh
      tr_surbox(2,23,block_no,cpu) = tr_neigh(2,6,block_no,cpu)
      !print *, "FACES:"
      !print *, "  5:", tr_neigh(1,5,block_no,cpu), tr_neigh(2,5,block_no,cpu)
      !print *, " 11:", tr_neigh(1,3,block_no,cpu), tr_neigh(2,3,block_no,cpu)
      !print *, " 13:", tr_neigh(1,1,block_no,cpu), tr_neigh(2,1,block_no,cpu)
      !print *, " 15:", tr_neigh(1,2,block_no,cpu), tr_neigh(2,2,block_no,cpu)
      !print *, " 17:", tr_neigh(1,4,block_no,cpu), tr_neigh(2,4,block_no,cpu)
      !print *, " 23:", tr_neigh(1,6,block_no,cpu), tr_neigh(2,6,block_no,cpu)



      ! middles of edges: 02,04,06,08, 10,12,16,18, 20,22,24,26
      if ((tr_surbox(1,5,block_no,cpu) .gt. 0) .and. (tr_surbox(2,5,block_no,cpu) .ge. 0)) then

        if (gr_bhTreeLrefine(block_no, cpu) .eq. gr_bhTreeLrefine(tr_surbox(1,5,block_no,cpu), tr_surbox(2,5,block_no,cpu))) then
          tr_surbox(1, 2,block_no,cpu) = tr_neigh(1, 3, tr_surbox(1,5,block_no,cpu), tr_surbox(2,5,block_no,cpu))
          tr_surbox(2, 2,block_no,cpu) = tr_neigh(2, 3, tr_surbox(1,5,block_no,cpu), tr_surbox(2,5,block_no,cpu))
          tr_surbox(1, 4,block_no,cpu) = tr_neigh(1, 1, tr_surbox(1,5,block_no,cpu), tr_surbox(2,5,block_no,cpu))
          tr_surbox(2, 4,block_no,cpu) = tr_neigh(2, 1, tr_surbox(1,5,block_no,cpu), tr_surbox(2,5,block_no,cpu))
          tr_surbox(1, 6,block_no,cpu) = tr_neigh(1, 2, tr_surbox(1,5,block_no,cpu), tr_surbox(2,5,block_no,cpu))
          tr_surbox(2, 6,block_no,cpu) = tr_neigh(2, 2, tr_surbox(1,5,block_no,cpu), tr_surbox(2,5,block_no,cpu))
          tr_surbox(1, 8,block_no,cpu) = tr_neigh(1, 4, tr_surbox(1,5,block_no,cpu), tr_surbox(2,5,block_no,cpu))
          tr_surbox(2, 8,block_no,cpu) = tr_neigh(2, 4, tr_surbox(1,5,block_no,cpu), tr_surbox(2,5,block_no,cpu))
          !print *, "DER 5:"
          !print *, "  2:", tr_neigh(1, 3, tr_surbox(1,5,block_no,cpu), tr_surbox(2,5,block_no,cpu)), tr_neigh(2, 3, tr_surbox(1,5,block_no,cpu), tr_surbox(2,5,block_no,cpu))
          !print *, "  4:", tr_neigh(1, 1, tr_surbox(1,5,block_no,cpu), tr_surbox(2,5,block_no,cpu)), tr_neigh(2, 1, tr_surbox(1,5,block_no,cpu), tr_surbox(2,5,block_no,cpu))
          !print *, "  6:", tr_neigh(1, 2, tr_surbox(1,5,block_no,cpu), tr_surbox(2,5,block_no,cpu)), tr_neigh(2, 2, tr_surbox(1,5,block_no,cpu), tr_surbox(2,5,block_no,cpu))
          !print *, "  8:", tr_neigh(1, 4, tr_surbox(1,5,block_no,cpu), tr_surbox(2,5,block_no,cpu)), tr_neigh(2, 4, tr_surbox(1,5,block_no,cpu), tr_surbox(2,5,block_no,cpu))
        endif
      endif

      if ((tr_surbox(1,11,block_no,cpu) .gt. 0) .and. (tr_surbox(2,11,block_no,cpu) .ge. 0)) then
        if (gr_bhTreeLrefine(block_no, cpu) .eq. gr_bhTreeLrefine(tr_surbox(1,11,block_no,cpu), tr_surbox(2,11,block_no,cpu))) then
          tr_surbox(1,10,block_no,cpu) = tr_neigh(1, 1, tr_surbox(1,11,block_no,cpu), tr_surbox(2,11,block_no,cpu))
          tr_surbox(2,10,block_no,cpu) = tr_neigh(2, 1, tr_surbox(1,11,block_no,cpu), tr_surbox(2,11,block_no,cpu))
          tr_surbox(1,12,block_no,cpu) = tr_neigh(1, 2, tr_surbox(1,11,block_no,cpu), tr_surbox(2,11,block_no,cpu))
          tr_surbox(2,12,block_no,cpu) = tr_neigh(2, 2, tr_surbox(1,11,block_no,cpu), tr_surbox(2,11,block_no,cpu))
          !print *, "DER 11:"
          !print *, " 10:", tr_neigh(1, 1, tr_surbox(1,11,block_no,cpu), tr_surbox(2,11,block_no,cpu)), tr_neigh(2, 1, tr_surbox(1,11,block_no,cpu), tr_surbox(2,11,block_no,cpu))
          !print *, " 12:", tr_neigh(1, 2, tr_surbox(1,11,block_no,cpu), tr_surbox(2,11,block_no,cpu)), tr_neigh(2, 2, tr_surbox(1,11,block_no,cpu), tr_surbox(2,11,block_no,cpu))
        endif
      endif
      if ((tr_surbox(1,17,block_no,cpu) .gt. 0) .and. (tr_surbox(2,17,block_no,cpu) .ge. 0)) then
        if (gr_bhTreeLrefine(block_no, cpu) .eq. gr_bhTreeLrefine(tr_surbox(1,17,block_no,cpu), tr_surbox(2,17,block_no,cpu))) then
          tr_surbox(1,16,block_no,cpu) = tr_neigh(1, 1, tr_surbox(1,17,block_no,cpu), tr_surbox(2,17,block_no,cpu))
          tr_surbox(2,16,block_no,cpu) = tr_neigh(2, 1, tr_surbox(1,17,block_no,cpu), tr_surbox(2,17,block_no,cpu))
          tr_surbox(1,18,block_no,cpu) = tr_neigh(1, 2, tr_surbox(1,17,block_no,cpu), tr_surbox(2,17,block_no,cpu))
          tr_surbox(2,18,block_no,cpu) = tr_neigh(2, 2, tr_surbox(1,17,block_no,cpu), tr_surbox(2,17,block_no,cpu))
          !print *, "DER 17:"
          !print *, " 16:", tr_neigh(1, 1, tr_surbox(1,17,block_no,cpu), tr_surbox(2,17,block_no,cpu)), tr_neigh(2, 1, tr_surbox(1,17,block_no,cpu), tr_surbox(2,17,block_no,cpu))
          !print *, " 18:", tr_neigh(1, 2, tr_surbox(1,17,block_no,cpu), tr_surbox(2,17,block_no,cpu)), tr_neigh(2, 2, tr_surbox(1,17,block_no,cpu), tr_surbox(2,17,block_no,cpu))
        endif
      endif
 
      if ((tr_surbox(1,23,block_no,cpu) .gt. 0) .and. (tr_surbox(2,23,block_no,cpu) .ge. 0)) then
        if (gr_bhTreeLrefine(block_no, cpu) .eq. gr_bhTreeLrefine(tr_surbox(1,23,block_no,cpu), tr_surbox(2,23,block_no,cpu))) then
          tr_surbox(1,20,block_no,cpu) = tr_neigh(1, 3, tr_surbox(1,23,block_no,cpu), tr_surbox(2,23,block_no,cpu))
          tr_surbox(2,20,block_no,cpu) = tr_neigh(2, 3, tr_surbox(1,23,block_no,cpu), tr_surbox(2,23,block_no,cpu))
          tr_surbox(1,22,block_no,cpu) = tr_neigh(1, 1, tr_surbox(1,23,block_no,cpu), tr_surbox(2,23,block_no,cpu))
          tr_surbox(2,22,block_no,cpu) = tr_neigh(2, 1, tr_surbox(1,23,block_no,cpu), tr_surbox(2,23,block_no,cpu))
          tr_surbox(1,24,block_no,cpu) = tr_neigh(1, 2, tr_surbox(1,23,block_no,cpu), tr_surbox(2,23,block_no,cpu))
          tr_surbox(2,24,block_no,cpu) = tr_neigh(2, 2, tr_surbox(1,23,block_no,cpu), tr_surbox(2,23,block_no,cpu))
          tr_surbox(1,26,block_no,cpu) = tr_neigh(1, 4, tr_surbox(1,23,block_no,cpu), tr_surbox(2,23,block_no,cpu))
          tr_surbox(2,26,block_no,cpu) = tr_neigh(2, 4, tr_surbox(1,23,block_no,cpu), tr_surbox(2,23,block_no,cpu))
          !print *, "DER 23:"
          !print *, " 20:", tr_neigh(1, 3, tr_surbox(1,23,block_no,cpu), tr_surbox(2,23,block_no,cpu)), tr_neigh(2, 3, tr_surbox(1,23,block_no,cpu), tr_surbox(2,23,block_no,cpu))
          !print *, " 22:", tr_neigh(1, 1, tr_surbox(1,23,block_no,cpu), tr_surbox(2,23,block_no,cpu)), tr_neigh(2, 1, tr_surbox(1,23,block_no,cpu), tr_surbox(2,23,block_no,cpu))
          !print *, " 24:", tr_neigh(1, 2, tr_surbox(1,23,block_no,cpu), tr_surbox(2,23,block_no,cpu)), tr_neigh(2, 2, tr_surbox(1,23,block_no,cpu), tr_surbox(2,23,block_no,cpu))
          !print *, " 26:", tr_neigh(1, 4, tr_surbox(1,23,block_no,cpu), tr_surbox(2,23,block_no,cpu)), tr_neigh(2, 4, tr_surbox(1,23,block_no,cpu), tr_surbox(2,23,block_no,cpu))
        endif
      endif
 
      ! corners: 01,03,07,09, 19,21,25,27
      if ((tr_surbox(1,2,block_no,cpu) .gt. 0) .and. (tr_surbox(2,2,block_no,cpu) .ge. 0)) then
        if (gr_bhTreeLrefine(block_no, cpu) .eq. gr_bhTreeLrefine(tr_surbox(1,2,block_no,cpu), tr_surbox(2,2,block_no,cpu))) then
          tr_surbox(1, 1,block_no,cpu) = tr_neigh(1, 1, tr_surbox(1, 2,block_no,cpu), tr_surbox(2, 2,block_no,cpu))
          tr_surbox(2, 1,block_no,cpu) = tr_neigh(2, 1, tr_surbox(1, 2,block_no,cpu), tr_surbox(2, 2,block_no,cpu))
          tr_surbox(1, 3,block_no,cpu) = tr_neigh(1, 2, tr_surbox(1, 2,block_no,cpu), tr_surbox(2, 2,block_no,cpu))
          tr_surbox(2, 3,block_no,cpu) = tr_neigh(2, 2, tr_surbox(1, 2,block_no,cpu), tr_surbox(2, 2,block_no,cpu))
          !print *, "CORNERS FROM 2:"
          !print *, "  1:", tr_neigh(1, 1, tr_surbox(1, 2,block_no,cpu), tr_surbox(2, 2,block_no,cpu)), tr_neigh(2, 1, tr_surbox(1, 2,block_no,cpu), tr_surbox(2, 2,block_no,cpu))
          !print *, "  3:", tr_neigh(1, 2, tr_surbox(1, 2,block_no,cpu), tr_surbox(2, 2,block_no,cpu)), tr_neigh(2, 2, tr_surbox(1, 2,block_no,cpu), tr_surbox(2, 2,block_no,cpu))
        endif
      endif
      if ((tr_surbox(1,8,block_no,cpu) .gt. 0) .and. (tr_surbox(2,8,block_no,cpu) .ge. 0)) then
        if (gr_bhTreeLrefine(block_no, cpu) .eq. gr_bhTreeLrefine(tr_surbox(1,8,block_no,cpu), tr_surbox(2,8,block_no,cpu))) then
          tr_surbox(1, 7,block_no,cpu) = tr_neigh(1, 1, tr_surbox(1, 8,block_no,cpu), tr_surbox(2, 8,block_no,cpu))
          tr_surbox(2, 7,block_no,cpu) = tr_neigh(2, 1, tr_surbox(1, 8,block_no,cpu), tr_surbox(2, 8,block_no,cpu))
          tr_surbox(1, 9,block_no,cpu) = tr_neigh(1, 2, tr_surbox(1, 8,block_no,cpu), tr_surbox(2, 8,block_no,cpu))
          tr_surbox(2, 9,block_no,cpu) = tr_neigh(2, 2, tr_surbox(1, 8,block_no,cpu), tr_surbox(2, 8,block_no,cpu))
          !print *, "CORNERS FROM 8:"
          !print *, "  7:", tr_neigh(1, 1, tr_surbox(1, 8,block_no,cpu), tr_surbox(2, 8,block_no,cpu)), tr_neigh(2, 1, tr_surbox(1, 8,block_no,cpu), tr_surbox(2, 8,block_no,cpu))
          !print *, "  9:", tr_neigh(1, 2, tr_surbox(1, 8,block_no,cpu), tr_surbox(2, 8,block_no,cpu)), tr_neigh(2, 2, tr_surbox(1, 8,block_no,cpu), tr_surbox(2, 8,block_no,cpu))
        endif
      endif
 
      if ((tr_surbox(1,20,block_no,cpu) .gt. 0) .and. (tr_surbox(2,20,block_no,cpu) .ge. 0)) then
        if (gr_bhTreeLrefine(block_no, cpu) .eq. gr_bhTreeLrefine(tr_surbox(1,20,block_no,cpu), tr_surbox(2,20,block_no,cpu))) then
          tr_surbox(1,19,block_no,cpu) = tr_neigh(1, 1, tr_surbox(1,20,block_no,cpu), tr_surbox(2,20,block_no,cpu))
          tr_surbox(2,19,block_no,cpu) = tr_neigh(2, 1, tr_surbox(1,20,block_no,cpu), tr_surbox(2,20,block_no,cpu))
          tr_surbox(1,21,block_no,cpu) = tr_neigh(1, 2, tr_surbox(1,20,block_no,cpu), tr_surbox(2,20,block_no,cpu))
          tr_surbox(2,21,block_no,cpu) = tr_neigh(2, 2, tr_surbox(1,20,block_no,cpu), tr_surbox(2,20,block_no,cpu))
          !print *, "CORNERS FROM 20:"
          !print *, " 19:", tr_neigh(1, 1, tr_surbox(1, 20,block_no,cpu), tr_surbox(2, 20,block_no,cpu)), tr_neigh(2, 1, tr_surbox(1, 20,block_no,cpu), tr_surbox(2, 20,block_no,cpu))
          !print *, " 21:", tr_neigh(1, 2, tr_surbox(1, 20,block_no,cpu), tr_surbox(2, 20,block_no,cpu)), tr_neigh(2, 2, tr_surbox(1, 20,block_no,cpu), tr_surbox(2, 20,block_no,cpu))
        endif
      endif
      if ((tr_surbox(1,26,block_no,cpu) .gt. 0) .and. (tr_surbox(2,26,block_no,cpu) .ge. 0)) then
        if (gr_bhTreeLrefine(block_no, cpu) .eq. gr_bhTreeLrefine(tr_surbox(1,26,block_no,cpu), tr_surbox(2,26,block_no,cpu))) then
          tr_surbox(1,25,block_no,cpu) = tr_neigh(1, 1, tr_surbox(1,26,block_no,cpu), tr_surbox(2,26,block_no,cpu))
          tr_surbox(2,25,block_no,cpu) = tr_neigh(2, 1, tr_surbox(1,26,block_no,cpu), tr_surbox(2,26,block_no,cpu))
          tr_surbox(1,27,block_no,cpu) = tr_neigh(1, 2, tr_surbox(1,26,block_no,cpu), tr_surbox(2,26,block_no,cpu))
          tr_surbox(2,27,block_no,cpu) = tr_neigh(2, 2, tr_surbox(1,26,block_no,cpu), tr_surbox(2,26,block_no,cpu))
          !print *, "CORNERS FROM 26:"
          !print *, " 25:", tr_neigh(1, 1, tr_surbox(1, 26,block_no,cpu), tr_surbox(2, 26,block_no,cpu)), tr_neigh(2, 1, tr_surbox(1, 26,block_no,cpu), tr_surbox(2, 26,block_no,cpu))
          !print *, " 27:", tr_neigh(1, 2, tr_surbox(1, 26,block_no,cpu), tr_surbox(2, 26,block_no,cpu)), tr_neigh(2, 2, tr_surbox(1, 26,block_no,cpu), tr_surbox(2, 26,block_no,cpu))
        endif
      endif

      !do i=1,6
      !  print *, "SURBOX N: ", i, cpu, block_no, tr_neigh(1,i,block_no,cpu), tr_neigh(2,i,block_no,cpu)
      !enddo

      !do i=1,27
      !  print *, "SURBOX S: ", i, cpu, block_no, tr_surbox(1,i,block_no,cpu), tr_surbox(2,i,block_no,cpu)
      !enddo


    enddo
  enddo

  return
end subroutine gr_bhFindNeighbours


