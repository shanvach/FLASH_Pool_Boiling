

#include "Flash.h"
#include "constants.h"
#include "Superlu.h"

subroutine gr_slugetSubsetBlks(slu_subset_flag,blockList_set,blockCount_set)


  use superlu_common, only : slu_max_FullyRefinedLevel, &
                             dofs_block, nloc_dofs

  use Grid_data, only : gr_meshMe,gr_meshNumProcs,gr_meshComm

  use Grid_interface, only : Grid_getMaxCommonRefinement, &
                             Grid_getListOfBlocks

  use Driver_interface, ONLY : Driver_abortFlash

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none
  integer, intent(in)  :: slu_subset_flag
  integer, intent(out) :: blockList_set(MAXBLOCKS),blockCount_set

  integer, save :: mgrid_solveLevelKPD

  
  ! Local Variables

  ! Initialize vars:
  blockCount_set   = 0
  blockList_set(:) = 0

  call RuntimeParameters_get('mgrid_solveLevelKPD', mgrid_solveLevelKPD)

  ! Obtain the subset of Blocks on the Processor:
  select case (slu_subset_flag)

  case(MAXUNIF_LEVEL)

     ! Define a refinement level to do the solve: 
     !call Grid_getMaxCommonRefinement(gr_MeshComm, slu_max_FullyRefinedLevel)
     slu_max_FullyRefinedLevel = mgrid_solveLevelKPD

     ! Make a list of blocks that belong to the level:
     call Grid_getListOfBlocks(REFINEMENT,blockList_set,blockCount_set, &
                               refinementLevel=slu_max_FullyRefinedLevel)

     print*,"SUPERLU SOLVE LEVEL:",slu_max_FullyRefinedLevel,"BlockCount",blockCount_set

  case default

     call Driver_abortFlash("gr_slugetSubsetBlks : slu_subset_flag not defined.")

  end select

  ! Allocate Matrix arrays:
  nloc_dofs    = dofs_block*blockCount_set

  return
end subroutine gr_slugetSubsetBlks
