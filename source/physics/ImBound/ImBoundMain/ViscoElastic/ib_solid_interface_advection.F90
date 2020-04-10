!===============================================================================
!!
!! subroutine ib_solid_interface_advection
!!
!! recontructs level set function for solid-fluid interface using the advected X and Y grid 
!!
!! X and Y grid initially overlaps with background Cartesian grid x and y 

!! Definition of the level set for solid-fluid interface remains the same  
!! in terms of X and Y throughout the simulation, where X(x,y), Y(x,y)
!===============================================================================

#include "constants.h"
#include "Flash.h"

         subroutine ib_solid_interface_advection(sd,sX,sY,&
                                   ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz)   
        implicit none
        !include 'mpif.h'
        real, dimension(:,:,:), intent(inout) :: sd,sX,sY

        real, intent(in)    :: dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2

        integer :: i,j,k
        !end header

           !!!!! reconsruct sd based on sX and sY
           ! below is a sample function for level set. It needs to be changed for a specific problem
             !sd = 0.d0
              k = 1
              do j = jy1-NGUARD, jy2+NGUARD
                 do i = ix1-NGUARD, ix2+NGUARD
                  sd(i,j,k) = - 0.2 + sqrt((sX(i,j,k)-0.1)**2 + &
                                         (sY(i,j,k))**2)
                 end do
              end do
         !print *,"ib_solid_interface_advection"
         end subroutine ib_solid_interface_advection
