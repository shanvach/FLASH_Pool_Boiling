!===============================================================================
!!
!! subroutine ib_levelset_linearprojection 
!!
!! extrapolating level set gradient in normal direction using 1st order upwind 
!!
!! args inout:: s       - levelset funtion to be projected
!! args in   :: so      - levelset at time n-1
!!           :: sn      - directional derivate of levelset 
!!           :: u,v     - normal vectors nx,ny
!!           :: dx,dy   - spacing in x and y directions
!!           :: ix1,ix2 - low,high block x indeces
!!           :: jy1,jy2 - low,high block y indeces 
!===============================================================================
        subroutine ib_levelset_linearprojection(s,sn,u,v,ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz)
        implicit none
        !include 'mpif.h'
        real, dimension(:,:,:), intent(inout) :: s
        real, dimension(:,:,:), intent(in)    :: u,v
        real, dimension(:,:,:), intent(in)    :: sn
        real, intent(in)    :: dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2

        integer :: i,j,k

        so(:,:,:) = s(:,:,:)

        !solve advection equation to get advected level set s with well defined sn  
        k = 1
        do j = jy1,jy2
          do i = ix1,ix2
          !normal vectors
          ul = u(i-1,j,k) 
          ur = u(i,j,k)   
          vl = v(i,j-1,k) 
          vr = v(i,j,k) 
          ! use dx/2 as dt to advect level set
          if(u(i,j,k).ge.0.d0.and.v(i,j,k).ge.0.d0) then
            s(i,j,k) = so(i,j,k) - dx/2.d0*u(i,j,k)*(so(i,j,k) - so(i-1,j,k)) / dx &
                                 - dx/2.d0*v(i,j,k)*(so(i,j,k) - so(i,j-1,k)) / dy &
                                 + dx/2.d0*sn(i,j,k)
          end if
          if(u(i,j,k).ge.0.d0.and.v(i,j,k).le.0.d0) then
            s(i,j,k) = so(i,j,k) - dx/2.d0*u(i,j,k)*(so(i,j,k) - so(i-1,j,k)) / dx &
                                 - dx/2.d0*v(i,j,k)*(so(i,j+1,k) - so(i,j,k)) / dy &
                                 + dx/2.d0*sn(i,j,k)
          end if
          if(u(i,j,k).le.0.d0.and.v(i,j,k).ge.0.d0) then
            s(i,j,k) = so(i,j,k) - dx/2.d0*u(i,j,k)*(so(i+1,j,k) - so(i,j,k)) / dx &
                                 - dx/2.d0*v(i,j,k)*(so(i,j,k) - so(i,j-1,k)) / dy &
                                 + dx/2.d0*sn(i,j,k)
          end if
          if(u(i,j,k).le.0.d0.and.v(i,j,k).le.0.d0) then
            s(i,j,k) = so(i,j,k) - dx/2.d0*u(i,j,k)*(so(i+1,j,k) - so(i,j,k)) / dx &
                                 - dx/2.d0*v(i,j,k)*(so(i,j+1,k) - so(i,j,k)) / dy &
                                 + dx/2.d0*sn(i,j,k)
          end if

          end do
        end do

        end subroutine ib_levelset_linearprojection
