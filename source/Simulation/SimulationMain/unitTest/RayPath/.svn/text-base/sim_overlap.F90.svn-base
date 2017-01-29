
#include "constants.h"
#include "Flash.h"

subroutine sim_overlap(bbox1, bbox2, overlap)

  real,dimension(LOW:HIGH,MDIM), intent(IN) :: bbox1,bbox2
  logical, intent(OUT) :: overlap

  logical :: noOverlap, temp


  overlap = .true.
  noOverlap = .false.

  do i = 1,NDIM
     if(overlap) then
        if(bbox1(LOW,i)<bbox2(LOW,i)) then
           temp=(bbox1(HIGH,i)<bbox2(LOW,i))
        else
           temp=(bbox2(HIGH,i)<bbox1(LOW,i))
        end if
     else
        temp=.false.
     end if
     noOverlap=noOverlap.or.temp
     overlap=.not.noOverlap
  
     if(overlap) then
        if(bbox1(HIGH,i)>bbox2(HIGH,i)) then
           temp=(bbox1(LOW,i)>bbox2(HIGH,i))
        else
           temp=(bbox2(LOW,i)>bbox1(HIGH,i))
        end if
     else
        temp = .false.
     end if
     
     noOverlap= noOverlap.or.temp
     overlap=.noOverlap
  end do

  return
end subroutine sim_overlap
  
