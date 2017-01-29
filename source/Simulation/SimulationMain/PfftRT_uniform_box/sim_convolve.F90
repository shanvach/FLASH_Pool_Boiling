subroutine sim_convolve(data_arr, limits, width, globalsize)

#include "Flash.h"
#include "constants.h"
#include "Pfft.h"

  use gr_pfftData, ONLY : pfft_outLen

  implicit none

  real, dimension(:), intent(INOUT) :: data_arr
  real, intent(IN) :: width
  integer, dimension(LOW:HIGH, MDIM), intent(IN) :: limits
  integer, dimension(MDIM), intent(IN) :: globalSize

  integer :: jx, jy, jz, ix, iy, iz, nx, ny, nz, indx
  real :: diag, qx, qy, qz, lx, ly, lz, dummy
  
  !*************************
  ! Here we convolve with a Fourier-space normalized Cartesian window function.
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
        if ( iy .eq. 1 ) then
            qy = 1.0
        else
            dummy = PI * ly / globalsize(2)
            qy = sin( 2.0 * dummy * width) / &
                   ( 2.0 * width * sin( dummy ) )
        endif

        do jx = 1, nx
            ix = limits(LOW, 2) + jx - 1  ! ix is the global index.  lx is
            lx = ix - 1                   ! the wave number (which starts at zero)
            if ( ix .eq. 1 ) then
                qx = 1.0
            else
                dummy = PI * lx / globalsize(1)
                qx = sin( 2.0 * dummy * width) / &
                    ( 2.0 * width * sin( dummy ) )
            endif

            do jz = 1, nz
                iz = limits(LOW, 1) + jz - 1  ! iz is the global index.  lz is
                lz = iz - 1                   ! the wave number (which starts at zero)
                if ( jz .eq. 1 ) then
                    qz = 1.0
                else
                    dummy = PI * lz / globalsize(3)
                    qz = sin( 2.0 * dummy * width) / &
                        ( 2.0 * width * sin( dummy ) )
                endif
            
                indx = 2*jz-1 + (jx - 1) * 2 * nz + (jy - 1) * 2 * nz * nx
        
                diag = qx * qy * qz
                data_arr(indx) = data_arr(indx) * diag
                data_arr(indx+1) = data_arr(indx+1) * diag
	
            end do
        end do
  end do
  
end subroutine sim_convolve
