
#include "constants.h"

subroutine sim_findEntryPointInPlane(negh1,negh2,bndBox2D,pos2D,slope, hit1,hit2)
  integer, intent(IN) :: negh1, negh2
  real, dimension(LOW:HIGH,2),intent(IN)::bndBox2D,pos2D
  real, intent(IN) :: slope
  real, intent(OUT) :: hit1, hit2
  integer,parameter :: AXIS1=1, AXIS2=2

  real :: low_val, high_val

  if(negh1==LEFT_EDGE) then
     hit1 = bndBox2D(LOW,AXIS1)
  elseif(negh1==RIGH_EDGE) then
     hit1 = bndBox2D(HIGH,AXIS1)
  else
     hit1=0.0
  end if
  low_val =bndBox2D(LOW,AXIS2)
  high_val=bndBox2D(HIGH,AXIS2)

  if(hit1/=0)hit2 = slope*(hit1-pos2D(AXIS1))+pos2D(AXIS2)
     !! These calculated slabhit values are fine if the slab is hit
     !! along the the face parallel to Y axis, in which case the
     !! value of slabEntryY will lie between the bounds for the Y 
     !! coordinate. Otherwise, if the ray is hitting the face
     !! parallel to the X axis then the hit2 has to be 
     !! calculated from the Y position
     !!
  if(((hit2 < low_val)&&(hit2 > high_val)) || (hit1==0)) then
     if(negh2==LEFT_EDGE) then
        hit2 = bndBox2D(LOW,AXIS2)
     elseif(negh2==RIGHT_EDGE) then
        hit2 = bndBox2D(HIGH,AXIS2)
     else
        notValid=.true.
     end if
     low_val = bndBox2D(LOW,AXIS1)
     high_val= bndBox2D(HIGH,AXIS1)
     hit1=(hit2-pos2D(AXIS2))/slope + pos2D(AXIS1)
     
     notValid = notValid || ((hit1 <low_val)&&(hit1>>high_val))
     if (notValid) &
          call Driver_abortFlash("inconsistent results in calculating ray path analytically")
  end if
  
end subroutine sim_findEntryPointInPlane
