subroutine Plasma_elPotential(ppold, distfunc, pp, dqn, dx, dy,ix1, ix2, jy1, jy2)

   use Plasma_data, only: pls_l2target

   implicit none
 
   real, dimension(:,:,:), intent(in) :: dqn, distfunc   
   real, dimension(:,:,:), intent(inout) :: pp, ppold  
   real, intent(in) :: dx, dy
   integer, intent(in) :: ix1, ix2, jy1, jy2

   real :: iter_diff, denominator 
   integer :: i,j
   !
   !iterations = 0
   iter_diff = pls_l2target + 1.0
   denominator = 0.0
   
   do while(iter_diff > pls_l2target)

     do j=jy1,jy2
        do i=ix1,ix2
           ppold(i,j,1) = pp(i,j,1)
        end do
     end do
     
     iter_diff = 0.0
     denominator = 0.0

     !all interior points
     do j=jy1+1,jy2-1
        do i=ix1+1,ix2-1
           if( distfunc(i,j,1) .ge. 0.0) then
           pp(i,j,1) = (0.5/((dx*dx)  + (dy*dy)))*                   &
                       ( (pp(i+1,j,1) + pp(i-1,j,1))*(dy*dy)       + &
                         (pp(i,j+1,1) + pp(i,j-1,1))*(dx*dx)       + &
                         (dqn(i,j,1))*(dx*dx)*(dy*dy)) 
           else
           pp(i,j,1) = (0.5/((dx*dx)  + (dy*dy)))*                   &
                       ( (pp(i+1,j,1) + pp(i-1,j,1))*(dy*dy)       + &
                         (pp(i,j+1,1) + pp(i,j-1,1))*(dx*dx))
           end if
        end do
     end do

     !homogeneous neumann boundary conditions
     !do j=jy1+1,jy2-1
        !left boundary
        !pp(ix1,j,1) = (1.0/(2*(dx*dx)  + (dy*dy)))*&
        !              ((pp(ix1+1,j,1)*(dy**2) + (pp(ix1,j+1,1) + pp(ix1,j-1,1))*(dx**2)) +&
        !              (dqn(ix1,j,1)/pls_epsilon0)*(dx*dx)*(dy*dy))
        !right boundary
        !pp(ix2,j,1) = (1.0/(2*(dx*dx)  + (dy*dy)))*&
        !              ((pp(ix2-1,j,1)*(dy**2) + (pp(ix2,j+1,1) + pp(ix2,j-1,1))*(dx**2)) +&
        !              (dqn(ix2,j,1)/pls_epsilon0)*(dx*dx)*(dy*dy))
     !end do

     !do i=ix1+1,ix2-1
        !lower boundary
        !pp(i,jy1,1) = (1.0/((dx*dx)  + 2*(dy*dy)))*&
        !              ((pp(i+1,jy1,1) + pp(i-1,jy1,1))*(dy**2) +&
        !               (dx**2)*pp(i,jy1+1,1) +&
        !               (dqn(i,jy1,1)/pls_epsilon0)*(dx*dx)*(dy*dy))
        !upper boundary
        !pp(i,jy2,1) = (1.0/((dx*dx)  + 2*(dy*dy)))*&
        !              (((pp(i+1,jy2,1) + pp(i-1,jy2,1))*(dy**2) + pp(i,jy2-1,1)*(dx**2)) +&
        !              (dqn(i,jy2,1)/pls_epsilon0)*(dx*dx)*(dy*dy))
     !end do

     !corners
     !pp(ix1,jy1,1) = (1.0/((dx**2 + dy**2)))*&
     !                ( pp(ix1+1,jy1,1)*(dy**2) + pp(ix1,jy1+1,1)*(dx**2) +&
     !                  (dqn(ix1,jy1,1)/pls_epsilon0)*(dx*dx)*(dy*dy)  )
                     !((pp(ix1+1,jy1,1)*(dy**2) + (pp(ix1,jy1+1,1))*(dx**2)) +&
                     ! (dqn(ix1,jy1,1)/pls_epsilon0)*(dx*dx)*(dy*dy))
     !pp(ix2,jy1,1) = (1.0/((dx**2 + dy**2)))*&
     !                ( pp(ix2-1,jy1,1)*(dy**2) + pp(ix2,jy1+1,1)*(dx**2) +&  
     !                  (dqn(ix2,jy1,1)/pls_epsilon0)*(dx*dx)*(dy*dy)  )
                     !((pp(ix2-1,jy1,1)*(dy**2) + (pp(ix2,jy1+1,1))*(dx**2)) +&
                     ! (dqn(ix2,jy1,1)/pls_epsilon0)*(dx*dx)*(dy*dy))
     !pp(ix1,jy2,1) = (1.0/((dx**2 + dy**2)))*&
     !                ( pp(ix1+1,jy2,1)*(dy**2) + pp(ix1,jy2-1,1)*(dx**2) +&
     !                  (dqn(ix1,jy2,1)/pls_epsilon0)*(dx*dx)*(dy*dy)  )
                     !((pp(ix1+1,jy2,1)*(dy**2) + (pp(ix1,jy2-1,1) )*(dx**2)) +&
                     !(dqn(ix1,jy2,1)/pls_epsilon0)*(dx*dx)*(dy*dy))
     !pp(ix2,jy2,1) = (1.0/((dx**2 + dy**2)))*&
     !                ( pp(ix2-1,jy2,1)*(dy**2) + pp(ix2,jy2-1,1)*(dx**2) +&
     !                  (dqn(ix2,jy2,1)/pls_epsilon0)*(dx*dx)*(dy*dy)  ) 
                     !((pp(ix2-1,jy2,1)*(dy**2) + (pp(ix2,jy2-1,1))*(dx**2)) +&
                     !(dqn(ix2,jy2,1)/pls_epsilon0)*(dx*dx)*(dy*dy))

     !compute L2norm of error
     do j=jy1+1,jy2-1
        do i=ix1+1,ix2-1
           iter_diff   = iter_diff   + (pp(i,j,1) - ppold(i,j,1))**2
           denominator = denominator + (ppold(i,j,1)*ppold(i,j,1)) 
        end do
     end do

     iter_diff = (iter_diff/denominator)
     iter_diff = iter_diff**0.5
     
   end do

end subroutine Plasma_elPotential
