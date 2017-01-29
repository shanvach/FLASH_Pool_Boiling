!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhBuildTreeBlock
!!
!! NAME
!!
!!  gr_bhBuildTree
!!
!!
!! SYNOPSIS
!!
!!   gr_bhBuildTreeBlock(
!!        integer,intent(in) :: block
!!        )
!!
!! DESCRIPTION
!!
!!   Build a tree for a specific block.
!!
!! ARGUMENTS
!!
!!  block    - ID of a block where the block-tree is constructed
!!
!!***


subroutine gr_bhBuildTreeBlock(idensvar, block)

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getCellCoords, &
    Grid_getBlkPtr, Grid_getDeltas
  use gr_bhInterface, ONLY : gr_bhGetTreeSize
  use gr_bhLocalInterface, ONLY : gr_bhGetTreePos
  use gr_bhData, ONLY: myPE => gr_bhTreeMyPE, gr_bhTreeLevels, gr_bhTreeArray, &
       gr_bhLocCoords, gr_bhTreeBS, gr_bhTreeParentTree
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer,intent(in) :: idensvar
  integer,intent(in) :: block
  integer       :: i, j, k, l, pos, fac, istat, level
  integer, dimension(2,MDIM)   :: blkLimits,blkLimitsGC
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  integer, DIMENSION(3) :: point
  integer       :: multi(1:gr_bhTreeLevels), mm
  real          :: del(MDIM)
  real          :: dvol, mass
  logical       :: gcell = .false.
  real          :: m = 0, xmc = 0, ymc = 0, zmc = 0
  



  ! allocate the tree array
  nullify(gr_bhTreeArray(myPE,block)%p)
  allocate(gr_bhTreeArray(myPE,block)%p(1:gr_bhGetTreeSize(gr_bhTreeLevels)), stat=istat)
  if (istat /= 0) then
    call Driver_abortFlash("could not allocate tree in gr_bhBuildTreeBlock.F90")
  endif

  ! fill tree with zeros
  do i = 1, gr_bhGetTreeSize(gr_bhTreeLevels)
    gr_bhTreeArray(myPE,block)%p(i) = 0
  enddo

  ! get information about coordinates of the block
  call Grid_getBlkIndexLimits(block,blkLimits,blkLimitsGC)
  call Grid_getCellCoords(IAXIS, block, CENTER, gcell, gr_bhLocCoords(:, IAXIS, block), gr_bhTreeBS)
  call Grid_getCellCoords(JAXIS, block, CENTER, gcell, gr_bhLocCoords(:, JAXIS, block), gr_bhTreeBS)
  call Grid_getCellCoords(KAXIS, block, CENTER, gcell, gr_bhLocCoords(:, KAXIS, block), gr_bhTreeBS)
  
  ! create the lowest level (and prepare the second lowest)
  call Grid_getBlkPtr(block,solnData,CENTER)
  call Grid_getDeltas(block,del)
  do k = 1,gr_bhTreeBS
    point(3) = k
    do j = 1,gr_bhTreeBS
      point(2) = j
      do i = 1,gr_bhTreeBS
        point(1) = i

        ! determine the multiindex
        do l = 1, gr_bhTreeLevels
          fac = 2**(gr_bhTreeLevels - l)
          multi(l) = 1 + mod((i-1),2*fac)/fac  &
                   + 2 *(mod((j-1),2*fac)/fac) &
                   + 4 *(mod((k-1),2*fac)/fac)
        enddo

        ! the lowest level - just write masses into the tree
        pos = gr_bhGetTreePos(gr_bhTreeLevels, multi)
        call Grid_getSingleCellVol(block, INTERIOR, point, dvol)
        mass = solnData(idensvar, i+blkLimits(LOW,IAXIS)-1, j+blkLimits(LOW,JAXIS)-1, k+blkLimits(LOW,KAXIS)-1)*dvol
        gr_bhTreeArray(myPE,block)%p(pos) = mass

        ! prepare the second lowest level (sum mass and m*(x,y,z) to the tree array)
        multi(gr_bhTreeLevels) = 0
        pos = gr_bhGetTreePos(gr_bhTreeLevels-1, multi)
        gr_bhTreeArray(myPE,block)%p(pos) = gr_bhTreeArray(myPE,block)%p(pos) + mass
        gr_bhTreeArray(myPE,block)%p(pos+1) = gr_bhTreeArray(myPE,block)%p(pos+1) + mass * gr_bhLocCoords(i, IAXIS, block)
        gr_bhTreeArray(myPE,block)%p(pos+2) = gr_bhTreeArray(myPE,block)%p(pos+2) + mass * gr_bhLocCoords(j, JAXIS, block)
        gr_bhTreeArray(myPE,block)%p(pos+3) = gr_bhTreeArray(myPE,block)%p(pos+3) + mass * gr_bhLocCoords(k, KAXIS, block)
        
      enddo
    enddo
  enddo

  ! second lowest level - normalization of the mass centers by the total mass is necessary
  do l = 1,gr_bhTreeLevels-1
    multi(l) = 1
  enddo
  multi(gr_bhTreeLevels) = 0
  do
    ! normalize mass center by the total mass in the branch
    pos = gr_bhGetTreePos(gr_bhTreeLevels-1, multi)
    gr_bhTreeArray(myPE,block)%p(pos+1) = gr_bhTreeArray(myPE,block)%p(pos+1) / gr_bhTreeArray(myPE,block)%p(pos)
    gr_bhTreeArray(myPE,block)%p(pos+2) = gr_bhTreeArray(myPE,block)%p(pos+2) / gr_bhTreeArray(myPE,block)%p(pos)
    gr_bhTreeArray(myPE,block)%p(pos+3) = gr_bhTreeArray(myPE,block)%p(pos+3) / gr_bhTreeArray(myPE,block)%p(pos)

    
    multi(gr_bhTreeLevels-1) = multi(gr_bhTreeLevels-1) + 1 ! update multi-index counter
    ! rotate multi-index counter
    do l = gr_bhTreeLevels-1,2,-1
      if (multi(l) > 8) then
        multi(l) = 1
        if (l > 1) multi(l-1) = multi(l-1) + 1
      endif
    enddo
    if (multi(1) > 8) exit 
  enddo
  

  ! create higher levels
  do level = gr_bhTreeLevels-1,1,-1
    
    ! set the initial multi-index
    do l = 1,gr_bhTreeLevels
      if (l <= level) then
        multi(l) = 1
      else
        multi(l) = 0
      endif
    enddo
    
    ! set accumulators to zero
    m = 0
    xmc = 0
    ymc = 0
    zmc = 0
    do
      ! sum mass and mass centres from the octet of branches 
      pos = gr_bhGetTreePos(level, multi)
      m = m + gr_bhTreeArray(myPE,block)%p(pos)
      xmc = xmc + gr_bhTreeArray(myPE,block)%p(pos+1) * gr_bhTreeArray(myPE,block)%p(pos)
      ymc = ymc + gr_bhTreeArray(myPE,block)%p(pos+2) * gr_bhTreeArray(myPE,block)%p(pos)
      zmc = zmc + gr_bhTreeArray(myPE,block)%p(pos+3) * gr_bhTreeArray(myPE,block)%p(pos)

      ! update multi-index counter
      multi(level) = multi(level) + 1

      if (multi(level) > 8) then
        ! write sum to the parent node
        mm = multi(level) ! backup multi(level)
        multi(level) = 0
        pos = gr_bhGetTreePos(level-1, multi)
        multi(level) = mm
        gr_bhTreeArray(myPE,block)%p(pos) = m
        gr_bhTreeArray(myPE,block)%p(pos+1) = xmc / m
        gr_bhTreeArray(myPE,block)%p(pos+2) = ymc / m
        gr_bhTreeArray(myPE,block)%p(pos+3) = zmc / m

        ! set accumulators to zero
        m = 0
        xmc = 0
        ymc = 0
        zmc = 0

      endif

      ! rotate the multi-index counter
      do l = level,2,-1
        if (multi(l) > 8) then
          multi(l) = 1
          if (l > 1) multi(l-1) = multi(l-1) + 1
        endif
      enddo
      if (multi(1) > 8) exit 
     
    enddo
  enddo

  gr_bhTreeParentTree(1:4, block,myPE) = gr_bhTreeArray(myPE,block)%p(1:4)

  return
end subroutine gr_bhBuildTreeBlock

