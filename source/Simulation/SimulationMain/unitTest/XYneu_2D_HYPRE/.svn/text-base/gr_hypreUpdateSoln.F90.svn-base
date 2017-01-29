!!****if* source/Simulation/SimulationMain/unitTest/XYneu_2D_HYPRE/gr_hypreUpdateSoln
!!
!!  NAME 
!!
!!  gr_hypreUpdateSoln
!!
!!  SYNOPSIS
!!
!!  call gr_hypreUpdateSoln (integer,intent(IN) :: iVar,
!!                           integer,intent(IN) :: blockCount,
!!                           integer,intent(IN) :: blockList (blockCount))
!!
!!  DESCRIPTION 
!!      This routine updates solution after solve (diffusion operation ,AX=B).
!!
!! ARGUMENTS
!!   iVar       : Variable on which the diffusion operatoion is performed (e.g TEMP_VAR)
!!   blockCount : The number of blocks in the list.   
!!   blockList  : The list of blocks on which the solution must be updated.  
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!!
!!   Uses HYPRE library.
!!   Solution is floored to 1.0e-16.
!!
!!***

!!REORDER(4): solnVec


subroutine gr_hypreUpdateSoln (iVar, blockCount, blockList)
  
  use Grid_interface,   ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr,      &       
       Grid_getBlkIndexLimits   
  use gr_hypreData,     ONLY : gr_hypreVecX, gr_hypreLower, &
       gr_hypreRefineMIN, gr_hypreUpper, gr_hypreVecB, gr_hypreFloor
  use Timers_interface, ONLY : Timers_start, Timers_stop    
  
  implicit none
  
#include "Flash.h"  
#include "constants.h"
#include "HYPREf.h"  
  
  integer,intent(IN) :: iVar
  integer,intent(IN) :: blockCount
  integer,intent(IN) :: blockList (blockCount)
  
  !! LOCAL VARIABLES
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  integer, dimension(2,MDIM):: blkLimitsGC, blkLimits  
  real :: values(1)
  integer :: blockID,part,level,var
  integer :: i, j, k, lb, ierr, pos(NDIM)
  integer :: ii
  real, allocatable, dimension(:) :: BoxVal
  integer :: datasize(MDIM)

    
  
  call Timers_start("gr_hypreUpdateSoln")    
  
  !!-----------------------------------------------------------------------
  !!     DEV: Update solution: small negative numbers are floored to 1.0e-16.
  !!-----------------------------------------------------------------------  
  
  call HYPRE_SStructVectorGather(gr_hypreVecX, ierr)  
  
  var = 0
  
  do lb = 1, blockCount     
     blockID = blockList(lb)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)    
     call Grid_getBlkPtr(blockID, solnVec)  
     call Grid_getBlkRefineLevel(blockID,level)                   
     
     part = level - gr_hypreRefineMIN
     
     datasize(1:MDIM)=blkLimits(HIGH,1:MDIM)-blkLimits(LOW,1:MDIM)+1
     
     allocate(BoxVal(product(dataSize(1:NDIM))))
     
     BoxVal = 0.0
     
     !! Use GetBoxValues more efficient then GetValues.
     call HYPRE_SStructVectorGetBoxValues(gr_hypreVecX, part,gr_hypreLower(lb,1:NDIM), &
          gr_hypreUpper(lb,1:NDIM), var, BoxVal(:), ierr)          
          
     do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)      
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
           do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)                      
              
              ii = (k-blkLimits(LOW,KAXIS)+1)                                +  &
                   (j-blkLimits(LOW,JAXIS))*dataSize(KAXIS)                  +  &
                   (i-blkLimits(LOW,IAXIS))*dataSize(KAXIS)*dataSize(JAXIS)
              
              solnVec(iVar,i,j,k) = BoxVal(ii)
                            
           end do
        end do
     end do

     BoxVal = 0.0
     call HYPRE_SStructVectorSetBoxValues(gr_hypreVecB, part,gr_hypreLower(lb,1:NDIM), &
          gr_hypreUpper(lb,1:NDIM), var, BoxVal(:), ierr) 
     
     deallocate (BoxVal)
     
     call Grid_releaseBlkPtr(blockID, solnVec)
  end do
  
  call Timers_stop("gr_hypreUpdateSoln") 
  
  return
  
end subroutine gr_hypreUpdateSoln
