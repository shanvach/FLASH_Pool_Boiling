Module superlu_common

  implicit none

  ! Local Dofs:
  integer, save :: nloc_dofs = 0

  ! Local and Total rows and columns:
  integer, save :: m_loc, m_glob, n_glob, fst_row

  ! Nonzeros:
  integer, parameter :: dofs_block = NXB*NYB*NZB
  integer, save :: nnz_allocate    
  integer, save :: nnz_loc, nnz_glob

  ! Matrix allocatable arrays:
  integer, save, allocatable, dimension(:) :: colind, rowptr
  real, save, allocatable, dimension(:) :: nzval  

  ! Grid Changed flag:
  logical, save :: gr_sluGridChanged = .false.

  ! Ref Level covering the whole domain:
  integer, save :: slu_max_FullyRefinedLevel = -1

  ! Neighbor Processors:
  integer, save, allocatable, dimension(:) :: neighProcsList
  integer, save :: neighProcsCount

  ! Global numeration startin index (0 based) for subset blks in the processor:
  integer, save :: blockCount_start_idx

  ! List of blockCount_start_idx and block_count_set of MyPE and neighbor Procs:
  integer, save, allocatable, dimension(:,:) :: neighProc_blkcnt
  
  ! List of subset blocks in local numeration of MyPE and neighbor Procs:
  integer, save, allocatable, dimension(:,:) :: neighProc_blkList

  ! List of global subset block numbers of MyPE and neighbor Procs:
  integer, save, allocatable, dimension(:,:) :: neighProc_blkList_idx


end Module superlu_common
