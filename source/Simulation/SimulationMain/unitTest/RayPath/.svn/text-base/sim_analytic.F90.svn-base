
#include "RayTrace.h"
#include "Flash.h"
#include "constants.h"

subroutine sim_analytic()
  use Simulation_data, ONLY : sim_numRay, sim_rayIncidence,sim_rayFinal,&
       sim_slabBndBox, sim_refract, sim_refractType, sim_domainBndBox
  use Grid_interface, ONLY : Grid_outsideBoundBox
  integer :: i,j

  real, dimension(MDIM) :: ,current_pos, next_pos
  real :: low_slope, high_slope
  integer, dimension(MDIM) :: negh
  logical :: outside
  integer, dimension(LEFT_EDGE:RIGHT_EDGE,LEFT_EDGE:RIGHT_EDGE,LOW:HIGH,2) :: &
       plane_endPoints
  
  logical :: hitSlab
  real, dimension(LOW,HIGH,2) :: bndBox2D, pos2D
  real :: slope1, slope2
  
  if(NDIM<2)call Driver_abortFlash("no refraction in a single dimension")

  plane_endPoints(LEFT_EDGE,LEFT_EDGE,LOW,1)  = HIGH
  plane_endPoints(LEFT_EDGE,LEFT_EDGE,LOW,2)  = LOW
  plane_endPoints(LEFT_EDGE,LEFT_EDGE,HIGH,1) = LOW
  plane_endPoints(LEFT_EDGE,LEFT_EDGE,HIGH,2) = HIGH

  plane_endPoints(LEFT_EDGE,CENTER,LOW,1)  = LOW
  plane_endPoints(LEFT_EDGE,CENTER,LOW,2)  = LOW
  plane_endPoints(LEFT_EDGE,CENTER,HIGH,1) = LOW
  plane_endPoints(LEFT_EDGE,CENTER,HIGH,2) = HIGH

  plane_endPoints(LEFT_EDGE,RIGHT_EDGE,LOW,1)  = LOW
  plane_endPoints(LEFT_EDGE,RIGHT_EDGE,LOW,2)  = LOW
  plane_endPoints(LEFT_EDGE,RIGHT_EDGE,HIGH,1) = HIGH
  plane_endPoints(LEFT_EDGE,RIGHT_EDGE,HIGH,2) = HIGH

  plane_endPoints(CENTER,LEFT_EDGE,LOW,1)  = HIGH
  plane_endPoints(CENTER,LEFT_EDGE,LOW,2)  = LOW
  plane_endPoints(CENTER,LEFT_EDGE,HIGH,1) = LOW
  plane_endPoints(CENTER,LEFT_EDGE,HIGH,2) = LOW

  plane_endPoints(CENTER,RIGHT_EDGE,LOW,1)  = LOW
  plane_endPoints(CENTER,RIGHT_EDGE,LOW,2)  = HIGH
  plane_endPoints(CENTER,RIGHT_EDGE,HIGH,1) = HIGH
  plane_endPoints(CENTER,RIGHT_EDGE,HIGH,2) = HIGH

  plane_endPoints(RIGHT_EDGE,LEFT_EDGE,LOW,1)  = LOW
  plane_endPoints(RIGHT_EDGE,LEFT_EDGE,LOW,2)  = LOW
  plane_endPoints(RIGHT_EDGE,LEFT_EDGE,HIGH,1) = HIGH
  plane_endPoints(RIGHT_EDGE,LEFT_EDGE,HIGH,2) = HIGH

  plane_endPoints(RIGHT_EDGE,CENTER,LOW,1)  = LOW
  plane_endPoints(RIGHT_EDGE,CENTER,LOW,2)  = HIGH
  plane_endPoints(RIGHT_EDGE,CENTER,HIGH,1) = HIGH
  plane_endPoints(RIGHT_EDGE,CENTER,HIGH,2) = HIGH

  plane_endPoints(RIGHT_EDGE,RIGHT_EDGE,LOW,1)  = LOW
  plane_endPoints(RIGHT_EDGE,RIGHT_EDGE,LOW,2)  = HIGH
  plane_endPoints(RIGHT_EDGE,RIGHT_EDGE,HIGH,1) = HIGH
  plane_endPoints(RIGHT_EDGE,RIGHT_EDGE,HIGH,2) = HIGH

  do i = 1,sim_numRay
     current_pos(1:MDIM) = sim_rayIncidence(RAY_IPOS:RAY_IPOS+MDIM-1,i)

     !! First determine which face is the ray going to hit
     call Grid_outsideBoundBox(current_pos,sim_slabBndBox,outside,negh)
     hitSlab=.true.

     !! First analyze the X-Y plane
     low_slope = sim_slabBndBox(plane_endPoints(neigh(IAXIS),&
                                neigh(JAXIS),LOW,2),JAXIS)-current_pos(JAXIS)/&
                                (plane_endPoints(neigh(IAXIS),neigh(JAXIS),LOW,1)-&
                                 current_pos(IAXIS))

     High_slope = sim_slabBndBox(plane_endPoints(neigh(IAXIS),&
                                neigh(JAXIS),HIGH,2),JAXIS)-current_pos(JAXIS)/&
                                (plane_endPoints(neigh(IAXIS),neigh(JAXIS),HIGH,1)-&
                                 current_pos(IAXIS))

     !! If only the slope of the ray in XY plane falls within the bounds is
     !! there any need to further process the ray, otherwise it is never going 
     !! to hit the slab

     if((sim_rayIncidence(RAY_XYSLOPE,i) < low_slope)|| &
          (sim_rayIncidence(RAY_XYSLOPE,i) > high_slope)) hitSlab=.false.

     if(hitSlab && (NDIM>2)) then

        !! If the dimensionality of the problem is 3, then the same
        !! test needs to be done for the YZ plane also
        
        low_slope = sim_slabBndBox(plane_endPoints(neigh(JAXIS),&
             neigh(KAXIS),LOW,2),KAXIS)-current_pos(KAXIS)/&
             (plane_endPoints(neigh(JAXIS),neigh(KAXIS),LOW,1)-&
             current_pos(JAXIS))
        
        High_slope = sim_slabBndBox(plane_endPoints(neigh(JAXIS),&
             neigh(KAXIS),HIGH,2),JAXIS)-current_pos(KAXIS)/&
             (plane_endPoints(neigh(JAXIS),neigh(KAXIS),HIGH,1)-&
             current_pos(JAXIS))
        
        !! If only the slope of the ray in XY plane falls within the bounds is
        !! there any need to further process the ray, otherwise it is never going 
        !! to hit the slab
         if((sim_rayIncidence(RAY_YZSLOPE,i) < low_slope).and.&
             (sim_rayIncidence(RAY_YZSLOPE,i) > high_slope))hitSlab=.false. 
     end if
     do j=1,2
        if(j==1) then
           slope1 = sim_rayIncidence(RAY_XYPLANE,i)
        else
           slope1 = sim_rayIncidence(RAY_YZPLANE,i)
        end if
        pos(1:2)=current_pos(j:j+1)

        if(hitSlab) then
           slope2=slope1*sim_refract
           bndBox2D(LOW:HIGH,1:2)=sim_slabBndBox(LOW:HIGH,j:j+1)
           call sim_findEntryPointInPlane(negh(j),negh(j+1),bndBox2D,&
                pos2D,slope1,next_pos(j),next_pos(j+1))
           
           pos(1:2)=next_pos(j:j+1)
           if(sim_refractType==RAY_LINEAR)then
              call sim_findLinearExitPointInPlane(bndBox2D,pos2D,slope2,&
                   next_pos(j),next_pos(j+1))
           else
              call sim_findCurvedExitPointInPlane(bndBox2D,pos2D,slope2,&
                   next_pos(j),next_pos(j+1))
           pos(1:2)=next_pos(j:j+1)
        end if
        bndBox2D(LOW:HIGH,1:2)=sim_domainBndBox(LOW:HIGH,j:j+1)
        call sim_findLinearExitPointInPlane(bndBox2D,pos2D,slope1,&
             next_pos(j),next_pos(j+1))
     end do
     sim_rayFinal(RAY_IPOS:RAY_IPOS+NDIM-1)=next_pos(1:NDIM)
  end do
end subroutine sim_analytic
