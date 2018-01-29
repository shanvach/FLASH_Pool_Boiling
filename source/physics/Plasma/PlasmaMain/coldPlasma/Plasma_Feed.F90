subroutine Plasma_Feed(constant_rate, noise, feed_rate, distfunc, ix1, ix2, jy1, jy2) 
   
   use Plasma_data, only: pls_pct_noise

   implicit none
  
   real, intent(in) :: constant_rate 
   real, dimension(:), intent(in) :: noise
   real, dimension(:,:,:), intent(in) :: distfunc 
 
   real, dimension(:,:,:), intent(inout) :: feed_rate

   integer, intent(in) :: ix1, ix2, jy1, jy2
   integer :: i,j

   real :: noise_shift

   call random_seed()

   do j=jy1,jy2
     do i=ix1,ix2
        !define boundary of plasma jet
        if( distfunc(i,j,1) .ge. 0.0) then
           !random number [0,1)
           call random_number(noise)
           !shift interval [-1,1)
           noise_shift = 2.0*noise(1) - 1.0
           !source rate for each node within the plasma boundary 
           !(noise based on 10% (arbitrary) value of constant rate 
           feed_rate(i,j,1) = constant_rate + pls_pct_noise*noise_shift*constant_rate 
        else
           feed_rate(i,j,1) = 0.0
        end if
      end do
   end do

end subroutine Plasma_Feed
