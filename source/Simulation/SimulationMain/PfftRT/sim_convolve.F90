subroutine sim_convolve(data_arr, kernel_arr, limits, globalsize)

#include "Flash.h"
#include "constants.h"
#include "Pfft.h"

  use gr_pfftData, ONLY : pfft_outLen
  use Simulation_data

  implicit none

  real, dimension(:), intent(INOUT) :: data_arr
  real, dimension(:), intent(IN) :: kernel_arr
  integer, dimension(LOW:HIGH, MDIM), intent(IN) :: limits
  integer, dimension(MDIM), intent(IN) :: globalSize

  integer :: jx, jy, jz, ix, iy, iz, nx, ny, nz, indx
  real :: fac, qx, qy, qz, lx, ly, lz, trigarg
  real, dimension(2) :: buf1, buf2

  !*************************
  ! Here we convolve with the kernel array whose FT is in kernel_arr
  !*************************

  !*************************
  ! Loop over local indices, compute Fourier-space Gaussian factor,
  ! multiply the corresponding element of data_arr.
  ! Keep in mind that the index order is now z,x,y, in consequence
  ! of the transpositions.
  !*************************

  nz = limits(HIGH, 1) - limits(LOW, 1) + 1
  nx = limits(HIGH, 2) - limits(LOW, 2) + 1
  ny = limits(HIGH, 3) - limits(LOW, 3) + 1
  do jy = 1, ny
    iy = limits(LOW, 3) + jy - 1  ! iy is the global index.  ly is
    ly = iy - 1                   ! the wave number (which starts at zero)

    do jx = 1, nx
      ix = limits(LOW, 2) + jx - 1  ! ix is the global index.  lx is
      lx = ix - 1                  ! the wave number (which starts at zero)

      do jz = 1, nz
        iz = limits(LOW, 1) + jz - 1  ! iz is the global index.  lz is
        lz = iz - 1                   ! the wave number (which starts at zero)
    
        indx = 2*jz-1 + (jx - 1) * 2 * nz + (jy - 1) * 2 * nz * nx
        trigarg = 2 * PI * ( lx * gauss_ctr_idx(1) / globalsize(1) + &
							 ly * gauss_ctr_idx(2) / globalsize(2) + &
							 lz * gauss_ctr_idx(3) / globalsize(3) )
        
        buf1(1) = cos(trigarg)
        buf1(2) = sin(trigarg)
        ! Shift kernel to origin
        buf2(1) = kernel_arr(indx) * buf1(1) - kernel_arr(indx+1) * buf1(2)
        buf2(2) = kernel_arr(indx) * buf1(2) + kernel_arr(indx+1) * buf1(1)
        ! Convolve shifted kernel with data
        buf1(1) = data_arr(indx) * buf2(1) - data_arr(indx+1) * buf2(2)
        buf1(2) = data_arr(indx) * buf2(2) + data_arr(indx+1) * buf2(1)

        data_arr(indx) = domain_volume * buf1(1)
        data_arr(indx+1) = domain_volume * buf1(2)
	
      end do
    end do
  end do
  
end subroutine sim_convolve
