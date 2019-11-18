!===============================================================================
!!
!! subroutine ib_levelset_constantprojection 
!!
!! extrapolating level set in direction normal to the interface using 1st order upwind
!!
!===============================================================================



      subroutine ib_levelset_constantprojection(s, sd, u, v, w)
        implicit none
        include 'mpif.h'
        integer i, j, k, npid, nop, nopx, nopy, nx, ny, nz, n
        integer nx1, nx2, ny1, ny2, nz1, nz2, dnx, dny, dnz
        double precision s, sd, sX, sY, u, v, w, x, y, z
        double precision xmx, xmn, ymx, ymn, zmn, zmx, dx, dy, dz
        double precision t, dt, ti, tf, cfl, delta_t
        dimension :: s(nx,ny,nz), u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
        dimension :: sd(nx,ny,nz), sX(nx,ny,nz), sY(nx,ny,nz)
        common/grid/xmx, xmn, ymx, ymn, zmn, zmx, dx, dy, dz, nx, ny, nz
        common/grid_index/nx1, nx2, ny1, ny2, nz1, nz2, dnx, dny, dnz
        common/mpi/npid, nop, nopx, nopy


        so = s
        
        do k = nz1, nz2

           
           do j = ny1-1, ny2+2


              do i = nx1-1, nx2+2


                 !cell face velocities
                 
                 ul = u(i-1,j,k) 
                 ur = u(i,j,k)   
                 vl = v(i,j-1,k) 
                 vr = v(i,j,k) 

                 
                 !!! use dx/2 as dt to advect level set
                 delta_t = dx/2.d0 

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

                 !!! use dx/2 as dt to advect level set


              end do
           end do
        end do
               
      end subroutine ib_levelset_constantprojection
