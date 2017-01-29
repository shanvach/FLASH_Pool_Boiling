!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhData
!!
!! NAME
!!
!!  gr_bhData
!!
!!
!! SYNOPSIS
!!
!!  use gr_bhData
!!
!!
!! DESCRIPTION
!!
!!  Variable declarations for the tree Poisson solver.
!!
!!   

module gr_bhData

  use Grid_interface, ONLY: GRID_PDE_BND_ISOLATED, GRID_PDE_BND_PERIODIC

  implicit none
#include "constants.h"

  integer, save              :: gr_bhTreeBS, gr_bhTreeLevels, gr_bhTreeZones
  integer, save              :: gr_bhTreeLoff(0:255)
  integer, save              :: gr_bhTreeMaxlnblocks
  integer, save              :: gr_bhTreeMyPE, gr_bhTreeNumProcs
  integer, save              :: gr_bhComm !replaces MPI_COMM_WORLD
  integer, save              :: gr_bhTreeLrefineMax
  real, save                 :: gr_bhTreeDcount(1:5)
  real, save                 :: gr_bhGravFac
  real, save                 :: gr_bhTreeLimangle, gr_bhTreeLimangle2
  !!DEV : The next two variables do not appear to be used or initialized anywhere. - KW
  real, save                 :: gr_bhTreeMincellmass, gr_bhTreeMaxcellmass
  real, save                 :: gr_bhLx, gr_bhLy, gr_bhLz
  real, save                 :: gr_bhEwaldXMax, gr_bhEwaldYMax, gr_bhEwaldZMax
  integer, save              :: gr_bhEwaldIsoFac
  character(len=MAX_STRING_LENGTH),save :: gr_bhEwaldFName
  logical, save              :: gr_bhEwaldAlwaysGenerate, gr_bhUseEwaldDecomp
  

  ! Supported boundary constants.
  integer, parameter :: GR_TREE_BND_ISOLATED = GRID_PDE_BND_ISOLATED !!  = 0
  integer, parameter :: GR_TREE_BND_PERIODIC = GRID_PDE_BND_PERIODIC !!  = 1

  ! Remember here the boundary condition for which Grid_solvePoisson is called.
  integer, save :: gr_bhBndType(6)

  ! array of pointers to trees 2D = (CPUs, blocks)
  type p_tree
    real, dimension(:), pointer :: p
  end type p_tree
  type(p_tree), save, dimension(:,:), allocatable :: gr_bhTreeArray

  ! array of pointers to messages 1D = (CPUs)
  type p_message
    real, dimension(:), pointer :: p
  end type p_message
  ! will be used like this in gr_bhExchangeTrees:
!!$  type(p_message), save, allocatable :: messages_send(:), messages_recv(:)

  ! levels of trees sent to different CPUs (block, toCPU)
  integer, save, allocatable :: gr_bhLocSentTreeLevels(:,:)
  ! levels of trees sent to different CPUs (block, fromCPU)
  integer, save, allocatable :: gr_bhLocRecvTreeLevels(:,:)

  ! Cell sizes (diagonals) (lrefine,MDIM)
  real, save, allocatable :: gr_bhTreeCellSize(:,:) ! is a function of lrefine and MDIM only
  ! Block sizes (diagonals) (lrefine)
  real, save, allocatable :: gr_bhTreeDiag2(:)

  ! Positions of all block centres (dim, block, cpu)
  real, save, allocatable :: gr_bhTreeBCen(:,:,:)

  ! Coordinates of all local blocks (gr_bhTreeBS+1, dim, block)
  real, save, allocatable :: gr_bhLocCoords(:,:,:)
  
  ! mass and mass centre positions of parent blocks (mass + MC position, block, cpu)
  ! filled also for LEAF blocks (in BuildTreeBlock) to communicate date for cases parent and child are at different CPU
  real, save, allocatable :: gr_bhTreeParentTree(:,:,:)
  real, save, allocatable :: gr_bhLocParentTree(:,:)

  ! amr tree structure - needed for parent blocks
  ! crated in treeComBlkProperties
  integer, save, allocatable :: gr_bhTreeNodetype(:,:)
  integer, save, allocatable :: gr_bhTreeLrefine(:,:)
  integer, save, allocatable :: gr_bhTreeChild(:,:,:,:)
  integer, save, allocatable :: gr_bhTreeNeigh(:,:,:,:)
  integer, save, allocatable :: gr_bhTreeLnblocks(:)

  ! stores indeces of leaf and parent blocks in a continous array (block_number, cpu)
  ! created in treeComBlkProperties
  integer, save, allocatable :: gr_bhTreeBlocklist(:, :)
  
  integer, save, allocatable :: gr_bhTreeFirstLevBlocks(:, :)
  integer, save :: gr_bhTreeNFirstLev

  ! correction table for Ewald summation - for periodic boundaries
  integer, save :: gr_bhEwaldFieldNx, gr_bhEwaldFieldNy, gr_bhEwaldFieldNz
  integer, save :: gr_bhEwaldSeriesN
  real, save, allocatable :: gr_bhTreeEwald(:,:,:)

  ! INTERACTION LISTS
  integer, save  :: gr_bhIlist

  ! 3^3 box of blocks surrounding a given block (cpu/tree, relative block position, block number, cpu)
  integer, save, allocatable :: gr_bhTreeSurbox(:,:,:,:)

  ! cell-node interaction list: list of indeces to the tree array
  integer, save, allocatable :: gr_bhTreeCn_ind(:,:,:)
  
  ! cell-cell interaction list:
  real, save, allocatable    :: gr_bhTreeCc_disti(:,:,:) ! list of inverted distances
  integer, save, allocatable :: gr_bhTreeCc_count(:,:,:) ! list of counts of the distances
  integer, save, allocatable :: gr_bhTreeCc_ind(:,:,:,:) ! appropriate indeces to the gr_bhTreeArray (masses)

  !! DEV: Next two unused - KW
  integer,save :: gr_bhCounter
  integer,save :: gr_bhTreeDebug_il

end module gr_bhData


