!!****if* source/physics/Diffuse/DiffuseMain/UG/diff_dotproduct
!!
!!  NAME 
!!
!!  diff_dotproduct
!!
!!  SYNOPSIS
!!
!!  call diff_dotproduct (A, B, dotprod, blkLimits, blkLimitsGC)
!!
!!  DESCRIPTION 
!!
!!
!! ARGUMENTS
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES:
!!  
!!
!!***

subroutine diff_dotproduct (vector1, vector2, dotprod, blkLimits, blkLimitsGC) 
  
  use Timers_interface,  ONLY : Timers_start, Timers_stop
  
  implicit none

 
#include "constants.h" 
  
  
  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC
  
  real, intent(IN)    :: vector1(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                                 blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                                 blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))  
  real, intent(IN)    :: vector2(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                                 blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                                 blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
  
  real, intent (OUT)  :: dotprod
  
  integer :: i, j, k
  
  dotprod = 0.  
  
  call Timers_start("dotproduct")
  
  do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
     do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
           dotprod = dotprod + vector1(i,j,k)*vector2(i,j,k)
        end do
     end do
  end do
  
  call Timers_stop("dotproduct")
  
end subroutine diff_dotproduct
