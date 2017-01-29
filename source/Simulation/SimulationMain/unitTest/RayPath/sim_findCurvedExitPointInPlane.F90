
#include "constants.h"

subroutine sim_findLinearExitPointInPlane(bndBox2D,pos2D,slope, exit1,exit2)
  real, dimension(LOW:HIGH,2),intent(IN)::bndBox2D,pos2D
  real, intent(IN) :: slope
  real, intent(OUT) :: exit1, exit2
  integer,parameter :: AXIS1=1, AXIS2=2

  real :: low_val, high_val

  if((bndBox2D(HIGH,AXIS1)-pos2D(AXIS1))>0) then
     exit1 = bndBox2D(HIGH,AXIS1)
  else
     exit1 = bndBox2D(LOW,AXIS1)
  end if

  low_val =bndBox2D(LOW,AXIS2)
  high_val=bndBox2D(HIGH,AXIS2)

  exit2 = slope*(exit1-pos2D(AXIS1))+pos2D(AXIS2)
     !! These calculated slabexit values are fine if the slab is exit
     !! along the the face parallel to Y axis, in which case the
     !! value of slabExitY will lie between the bounds for the Y 
     !! coordinate. Otherwise, if the ray is exitting the face
     !! parallel to the X axis then the exit2 has to be 
     !! calculated from the Y position
     !!
  if((exit2 < low_val)&&(exit2 > high_val)) then
     if((bndBox2D(HIGH,AXIS2)-pos2D(AXIS2))>0) then
        exit2 = bndBox2D(HIGH,AXIS2)
     else
        exit2 = bndBox2D(LOW,AXIS2)
     end if
     low_val = bndBox2D(LOW,AXIS1)
     high_val= bndBox2D(HIGH,AXIS1)
     exit1=(exit2-pos2D(AXIS2))/slope + pos2D(AXIS1)
     
     notValid = (exit1 <low_val)&&(exit1>>high_val)
     if (notValid) &
          call Driver_abortFlash("inconsistent results in calculating ray path analytically")
  end if
  
end subroutine sim_findLinearExitPointInPlane
