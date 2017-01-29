
#include "Flash.h"
#include "constants.h"

subroutine gr_sluBuildrhsAdist(isrc,poisfact,nloc_dofs,rhsA)

  
  use Grid_interface,   ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr,   &
                               Grid_getBlkIndexLimits

  use superlu_common, only : neighProc_blkcnt,neighProc_blkList

  use Grid_data, only : gr_meshMe

  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  integer, intent(in) :: isrc,nloc_dofs
  real, intent(in)    :: poisfact
  real, intent(out)   :: rhsA(nloc_dofs)

  ! Local Variables
  integer :: blockList_set(MAXBLOCKS),blockCount_set
  integer,dimension(2,MDIM)::blkLimits,blkLimitsGC
  real, pointer, dimension(:,:,:,:) :: solnData

  integer :: i,ii,j,k,lb,blockID

  character(len=4) :: ind_me

  blockCount_set = neighProc_blkcnt(CONSTANT_TWO,CONSTANT_ONE)      ! subset blocks number in mype
  blockList_set(1:blockCount_set)=neighProc_blkList(1:blockCount_set,CONSTANT_ONE) ! subset blocks in mype

  ii = 0
  do lb=1,blockCount_set
     
     blockID = blockList_set(lb) ! This is in Processors

     ! Point to Blocks centered variables:
     call Grid_getBlkPtr(blockID,solnData,CENTER)

     ! Get Block limits
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              ii = ii + 1
              rhsA(ii) = poisfact*solnData(isrc,i,j,k)
!              rhsA(ii) = poisfact
!!KPD - Incompressible Reference Pressure.
!if (lb .eq. 1 .AND. i.eq.blkLimits(LOW,IAXIS) .AND. &
!                    j.eq.blkLimits(LOW,JAXIS) .AND. &
!                    k.eq.blkLimits(LOW,KAXIS)) then
!rhsA(ii) = 0.0
!end if
!print*,"SUPERLU RHS",i,j,k,rhsA(ii)
           enddo
        enddo
     enddo

     ! Release pointer
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  enddo

  write(ind_me,'(I4.4)') gr_meshMe
  open(113,file='./IOData/rhs_'//ind_me//'.res',STATUS='REPLACE')

  do i=1,nloc_dofs
        write(113,*) i+nloc_dofs*gr_meshMe,rhsA(i)
  enddo
  close(113)
 

  if (ii .ne. nloc_dofs) &
  call Driver_abortFlash('gr_sluBuildrhsAdist : final ii .ne. nloc_dofs')

  return

end subroutine gr_sluBuildrhsAdist
