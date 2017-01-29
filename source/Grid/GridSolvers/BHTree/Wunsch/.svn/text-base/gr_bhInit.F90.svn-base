!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhInit
!!
!! NAME
!!
!!  gr_bhInit
!!
!! 
!! SYNOPSIS
!!
!!  call gr_bhInit()
!!
!!
!! DESCRIPTION
!!
!!  Initialize the tree Poisson solver.  Read in any of the
!!  runtime parameters for this solver.  All solver common data
!!  is stored in the gr_bhData module.
!!
!!
!!
!!***

subroutine gr_bhInit()

  use Grid_data, ONLY : gr_geometry, gr_isolatedBoundaries, &
       gr_meshComm, gr_meshMe, gr_meshNumProcs
  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use gr_bhData, ONLY : gr_bhTreeLevels, &
    gr_bhTreeArray, gr_bhTreeBS, gr_bhTreeLimAngle2, &
    gr_bhLocCoords, gr_bhTreeLoff, gr_bhTreeMyPE, &
    gr_bhTreeBCen, gr_bhTreeLrefine, gr_bhTreeCellSize, &
    gr_bhTreeDiag2, gr_bhLocRecvTreeLevels, &
    gr_bhLocSentTreeLevels, gr_bhTreeNodetype, &
    gr_bhComm, gr_bhTreeNumProcs, gr_bhTreeLnblocks, &
    gr_bhTreeParentTree, gr_bhLocParentTree, gr_bhTreeChild, &
    gr_bhTreeBlocklist, gr_bhTreeNeigh, gr_bhTreeSurBox, &
    gr_bhTreeCc_count, gr_bhTreeCc_disti, gr_bhTreeCn_ind, gr_bhTreeCc_ind, &
    gr_bhIlist, gr_bhGravFac, gr_bhTreeLrefineMax, gr_bhTreeLimAngle, &
    gr_bhTreeNFirstLev, gr_bhTreeFirstLevBlocks



  use tree, ONLY : nodetype, lrefine, child, mchild, maxblocks_tr, mfaces
  use Logfile_interface, ONLY : Logfile_stamp
  use gr_bhInterface, ONLY : gr_bhInitTemplates, gr_bhEwaldField
 
  implicit none
  
#include "Flash.h"
#include "Flash_mpi.h"
#include "constants.h"

  integer :: istat
  integer :: i, j, k, ni, nj, nk, hi, hj, hk
  integer :: nblockx, nblocky, nblockz
  integer :: max_cc_dist, max_cc_ind, max_cn_ind

  real :: mcs(MDIM)
  character(len=32) :: fname    !! DEV: unused
  character(len=MAX_STRING_LENGTH) :: strBuff


  gr_bhTreeMyPE = gr_meshMe
  gr_bhTreeNumProcs = gr_meshNumProcs
  gr_bhComm = gr_meshComm
  call RuntimeParameters_get("gr_bhTreeLimangle", gr_bhTreeLimangle)
  call RuntimeParameters_get("gr_bhIlist", gr_bhIlist)
  call RuntimeParameters_get("lrefine_max", gr_bhTreeLrefineMax)
  gr_bhTreeLimAngle2 = gr_bhTreeLimangle**2

  ! Check if we support the requested grid geometry.
  if ((NDIM /= 3) .or. (gr_geometry /= CARTESIAN)) then
     call Driver_abortFlash ('[gr_bhInit] ERROR: tree Poisson solver works only in 3D Cartesian geometry')
  endif

  ! check if NBX = NBY = NBZ
  if ((NXB /= NYB) .or. (NXB /= NZB) .or. (NYB /= NZB)) then
     call Driver_abortFlash ('[gr_bhInit] ERROR: must be NXB = NYB = NZB')
  endif
  
  ! set the tree parameters: gr_bhTreeBS, gr_bhTreeLevels
  gr_bhTreeBS = NXB
  gr_bhTreeLevels = 0 ! gr_bhTreeLevels = ln_2 (NBX)
  i = gr_bhTreeBS
  do
    if (mod(i, 2) == 0) then
      i = i / 2
      gr_bhTreeLevels = gr_bhTreeLevels + 1
    else 
      exit 
    endif
    if (i == 1) exit
  enddo
  do i = 0, gr_bhTreeLevels
    gr_bhTreeLoff(i) = 1 + 4*(8**i - 1)/7 ! offset of the level
  enddo

  ! cell size
  allocate(gr_bhTreeCellSize(1:gr_bhTreeLrefineMax,MDIM), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeCellSize in gr_bhInit.F90")
  call Grid_getMinCellSizes(mcs)
  gr_bhTreeCellSize(gr_bhTreeLrefineMax,:) = mcs
  do i = gr_bhTreeLrefineMax-1,1,-1
    gr_bhTreeCellSize(i,:) = 2*gr_bhTreeCellSize(i+1,:)
  enddo

  ! 3D diagonal of blocks
  allocate(gr_bhTreeDiag2(1:gr_bhTreeLrefineMax+gr_bhTreeLevels), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeDiag2 in gr_bhInit.F90")
  do i = 1, gr_bhTreeLrefineMax
    gr_bhTreeDiag2(i) = (NXB*gr_bhTreeCellSize(i,IAXIS))**2 &
                + (NYB*gr_bhTreeCellSize(i,JAXIS))**2 &
                + (NZB*gr_bhTreeCellSize(i,KAXIS))**2
  enddo
  do i = gr_bhTreeLrefineMax+1, gr_bhTreeLrefineMax+gr_bhTreeLevels
    gr_bhTreeDiag2(i) = gr_bhTreeDiag2(i-1)/4.0
  enddo

      
  allocate(gr_bhTreeArray(0:gr_bhTreeNumProcs-1,1:MAXBLOCKS), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeArray in gr_bhInit.F90")

  allocate(gr_bhLocSentTreeLevels(MAXBLOCKS, 0:gr_bhTreeNumProcs-1), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhLocSentTreeLevels in gr_bhInit.F90")
  allocate(gr_bhLocRecvTreeLevels(MAXBLOCKS, 0:gr_bhTreeNumProcs-1), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhLocRecvTreeLevels in gr_bhInit.F90")
  
  allocate(gr_bhLocCoords(gr_bhTreeBS+1,MDIM, MAXBLOCKS), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhLocCoords in gr_bhInit.F90")
  allocate(gr_bhTreeBCen(MDIM, MAXBLOCKS,0:gr_bhTreeNumProcs-1), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeBCen in gr_bhInit.F90")

  allocate(gr_bhTreeParentTree(4, MAXBLOCKS, 0:gr_bhTreeNumProcs-1), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeParentTrees in gr_bhInit.F90")
  allocate(gr_bhLocParentTree(4, MAXBLOCKS), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhLocParentTrees in gr_bhInit.F90")
  
  allocate(gr_bhTreeNodetype(maxblocks_tr, 0:gr_bhTreeNumProcs-1), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeNodetype in gr_bhInit.F90")
  allocate(gr_bhTreeLrefine(maxblocks_tr, 0:gr_bhTreeNumProcs-1), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeLrefine in gr_bhInit.F90")
  allocate(gr_bhTreeChild(2, mchild, maxblocks_tr, 0:gr_bhTreeNumProcs-1), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeChild in gr_bhInit.F90")
  allocate(gr_bhTreeLnblocks(0:gr_bhTreeNumProcs-1), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeLnblocks in gr_bhInit.F90")

  allocate(gr_bhTreeBlocklist(1:1, 0:gr_bhTreeNumProcs-1), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeBlockList in gr_bhInit.F90")

  call RuntimeParameters_get("nblockx", nblockx)
  call RuntimeParameters_get("nblocky", nblocky)
  call RuntimeParameters_get("nblockz", nblockz)
  gr_bhTreeNFirstLev = nblockx*nblocky*nblockz
  allocate(gr_bhTreeFirstLevBlocks(2, gr_bhTreeNFirstLev), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeFirstLevBlocks in gr_bhInit.F90")

  ! INTERACTION LISTS
  if (gr_bhIlist .ne. 0) then
    allocate(gr_bhTreeNeigh(2, mfaces, maxblocks_tr, 0:gr_bhTreeNumProcs-1), stat=istat)
    if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeNeigh in gr_bhInit.F90")

    allocate(gr_bhTreeSurbox(2, 27, MAXBLOCKS, 0:gr_bhTreeNumProcs-1), stat=istat)
    if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeSurbox in gr_bhInit.F90")

    call gr_bhInitTemplates(0, max_cc_dist, max_cc_ind, max_cn_ind)
    if (gr_bhTreeMyPE == MASTER_PE) then
       write (strBuff, '("max_cc_dist=", i6, " max_cc_inx=", i6, " max_cn_ind=", i6)') &
             max_cc_dist, max_cc_ind, max_cn_ind
       call Logfile_stamp(strBuff, "[BHTree]")
    end if

    ! cell-node interaction list: list of indeces to the tree array
    allocate(gr_bhTreeCn_ind(1:gr_bhTreeBS**3,1:max_cn_ind, 27), stat=istat)
    if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeCn_ind in gr_bhInit.F90")
    
    ! cell-cell interaction list:
    allocate(gr_bhTreeCc_disti(1:gr_bhTreeBS**3,1:max_cc_dist, 27), stat=istat)
    if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeCc in gr_bhInit.F90")
    allocate(gr_bhTreeCc_count(1:gr_bhTreeBS**3,1:max_cc_dist, 27), stat=istat)
    if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeCc_count in gr_bhInit.F90")
    allocate(gr_bhTreeCc_ind(1:gr_bhTreeBS**3,1:max_cc_dist,1:max_cc_ind, 27), stat=istat)
    if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhTreeCc_ind in gr_bhInit.F90")
 
    call gr_bhInitTemplates(1, max_cc_dist, max_cc_ind, max_cn_ind)

  endif

  ! PERIODIC BOUNDARIES
  call gr_bhEwaldField()
  
  ! to ensure all allocations are made before the code continues
  call MPI_Barrier(gr_bhComm, i) 
  
  return
end subroutine gr_bhInit


