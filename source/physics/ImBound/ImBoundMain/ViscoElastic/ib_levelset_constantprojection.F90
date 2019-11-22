!===============================================================================
!!
!! subroutine ib_levelset_constantprojection 
!!
!! extrapolating level set in direction normal to the interface using 1st order upwind
!!
!! args inout:: s       - levelset funtion to be projected
!! args in   :: so      - levelset function s at time n-1
!!           :: u,v     - normal vectors nx,ny
!!           :: dx,dy   - spacing in x and y directions
!!           :: ix1,ix2 - low,high block x indeces
!!           :: jy1,jy2 - low,high block y indeces 
!===============================================================================
        subroutine ib_levelset_constantprojection(s,so,u,v,dx,dy,ix1,ix2,jy1,jy2)
        implicit none
        real, dimension(:,:,:), intent(inout) :: s
        real, dimension(:,:,:), intent(in)    :: so,u,v
        real, intent(in)    :: dx,dy
        integer, intent(in) :: ix1,ix2,jy1,jy2

        integer :: i,j,k
        real    :: ul,ur,vl,vr,delta_t

        so(:,:,:) = s(:,:,:)
        
        k = 1

        !use dx/2 as dt to advect level set
        delta_t = dx/2.d0

        do j = jy1,jy2
          do i = ix1,ix2
          !normal vectors
          ul = u(i-1,j,k) 
          ur = u(i,j,k)   
          vl = v(i,j-1,k) 
          vr = v(i,j,k) 
          !use dx/2 as dt to advect level set
          if(u(i,j,k).ge.0.d0.and.v(i,j,k).ge.0.d0) then
            s(i,j,k) = so(i,j,k) - delta_t*u(i,j,k)*(so(i,j,k) - so(i-1,j,k)) / dx &
                                 - delta_t*v(i,j,k)*(so(i,j,k) - so(i,j-1,k)) / dy
          end if
          if(u(i,j,k).ge.0.d0.and.v(i,j,k).le.0.d0) then
            s(i,j,k) = so(i,j,k) - delta_t*u(i,j,k)*(so(i,j,k) - so(i-1,j,k)) / dx &
                                 - delta_t*v(i,j,k)*(so(i,j+1,k) - so(i,j,k)) / dy
          end if
          if(u(i,j,k).le.0.d0.and.v(i,j,k).ge.0.d0) then
            s(i,j,k) = so(i,j,k) - delta_t*u(i,j,k)*(so(i+1,j,k) - so(i,j,k)) / dx &
                                 - delta_t*v(i,j,k)*(so(i,j,k) - so(i,j-1,k)) / dy
          end if
          if(u(i,j,k).le.0.d0.and.v(i,j,k).le.0.d0) then
            s(i,j,k) = so(i,j,k) - delta_t*u(i,j,k)*(so(i+1,j,k) - so(i,j,k)) / dx &
                                 - delta_t*v(i,j,k)*(so(i,j+1,k) - so(i,j,k)) / dy
          end if
          end do
        end do
               
      end subroutine ib_levelset_constantprojection
