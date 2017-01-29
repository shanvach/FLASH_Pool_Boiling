Module slu_interface

#include "Flash.h"

  implicit none

  interface
     subroutine gr_slugetSubsetBlks(slu_subset_flag,blockList_set,blockCount_set)
       implicit none
       integer, intent(in)  :: slu_subset_flag
       integer, intent(out) :: blockList_set(MAXBLOCKS),blockCount_set
     end subroutine gr_slugetSubsetBlks
  end interface

  interface
     subroutine gr_slugetSubsetStartidx(blockCount_set) 
       implicit none
       integer, intent(in) :: blockCount_set
     end subroutine gr_slugetSubsetStartidx
  end interface

  interface
     subroutine gr_slugetSubsetNeighProcs(blockList_set,blockCount_set,diagblk_flag)
       implicit none
       integer, intent(in) :: blockList_set(MAXBLOCKS),blockCount_set
       integer, intent(in) :: diagblk_flag
     end subroutine gr_slugetSubsetNeighProcs
  end interface

  interface
     subroutine gr_slugetSubsetNeighLists(blockList_set,blockCount_set)
       implicit none
       integer, intent(in) :: blockList_set(MAXBLOCKS),blockCount_set
     end subroutine gr_slugetSubsetNeighLists
  end interface

  interface
     subroutine gr_sluBuildAdist(A, m, n)
       !use superlu_common, only : m_loc , nnz_allocate
       use superlu_mod
       implicit none
       integer(superlu_ptr), intent(in) :: A
       integer, intent(inout) :: m, n
       !integer, intent(inout) :: colind(nnz_allocate), rowptr(m_loc+1)  
       !real, intent(inout) :: nzval(nnz_allocate)
     end subroutine gr_sluBuildAdist
  end interface

  interface
     subroutine gr_sluBuildrhsAdist(isrc,poisfact,nloc_dofs,rhsA)
       implicit none
       integer, intent(in) :: isrc,nloc_dofs
       real, intent(in)    :: poisfact
       real, intent(out)   :: rhsA(nloc_dofs)
     end subroutine gr_sluBuildrhsAdist
  end interface

  interface
     subroutine gr_sluDumpSolndist(isoln,nloc_dofs,solA)
       implicit none
       integer, intent(in) :: isoln,nloc_dofs
       real, intent(in)   :: solA(nloc_dofs)
     end subroutine gr_sluDumpSolndist
  end interface

End Module slu_interface
