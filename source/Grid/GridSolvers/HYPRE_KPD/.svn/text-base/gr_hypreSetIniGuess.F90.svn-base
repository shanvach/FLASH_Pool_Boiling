!!****if* source/Grid/GridSolvers/HYPRE_KPD/gr_hypreSetIniGuess
!!
!!  NAME 
!!
!!  gr_hypreSetIniGuess
!!
!!  SYNOPSIS
!!
!!  call gr_hypreSetIniGuess (integer,intent(IN) :: iVar,
!!                            integer,intent(IN) :: blockCount,
!!                            integer,intent(IN) :: blockList (blockCount))
!!
!!  DESCRIPTION 
!!      This routine sets initial guess for AX=B.
!!
!! ARGUMENTS
!!   iVar       : Variable on which the diffusion operatorion is performed (e.g TEMP_VAR)
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
!!
!!***

!!REORDER(4): solnVec


subroutine gr_hypreSetIniGuess (iVar, blockCount, blockList)
  
  use gr_hypreData,     ONLY : gr_hypreVecX, gr_hypreLower, &
                               gr_hypreRefineMIN, gr_hypreUpper
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_getBlkIndexLimits, Grid_getBlkRefineLevel
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
  integer :: blockID,part,level,var
  integer :: i, j, k, lb, ierr, pos(NDIM)
  integer :: ii
  real, allocatable, dimension(:) :: BoxVal
  integer :: datasize(MDIM)
  
  call Timers_start("gr_hypreSetIniGuess")     
  
  var = 0
  
  do lb = 1, blockCount  
     blockID = blockList(lb)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)    
     call Grid_getBlkPtr(blockID,solnVec)  

     call Grid_getBlkRefineLevel(blockID,level)       
     part = level - gr_hypreRefineMIN
     
     datasize(1:MDIM)=blkLimits(HIGH,1:MDIM)-blkLimits(LOW,1:MDIM)+1     
     allocate(BoxVal(product(dataSize(1:NDIM))))     
     BoxVal = 0.0
     
     do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)      
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
           do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)                            
              
              ii = (k-blkLimits(LOW,KAXIS)+1)                                +  &
                   (j-blkLimits(LOW,JAXIS))*dataSize(KAXIS)                  +  &
                   (i-blkLimits(LOW,IAXIS))*dataSize(KAXIS)*dataSize(JAXIS)            
              
              BoxVal(ii) = solnVec(iVar,i,j,k)               
              
           end do
        end do
     end do
     
     !! Initial guess
     call HYPRE_SStructVectorSetBoxValues(gr_hypreVecX, part,gr_hypreLower(lb,1:NDIM), &
          gr_hypreUpper(lb,1:NDIM), var, BoxVal(:), ierr)               
     
     deallocate (BoxVal)
     call Grid_releaseBlkPtr(blockID, solnVec)    
     
  end do
  
  call HYPRE_SStructVectorAssemble(gr_hypreVecX, ierr)  
  
  
  call Timers_stop("gr_hypreSetIniGuess") 
  
  return
  
end subroutine gr_hypreSetIniGuess
