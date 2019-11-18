!===============================================================================
!!
!! subroutine ib_dynamic_grid_advection 
!!
!! extrapolating a scalar field from inside to the outside of a region defined by sd
!!
!===============================================================================



      subroutine ib_dynamic_grid_advection(sn, sd, stest, adfx, adfy, cpt)
        implicit none
        include 'mpif.h'

        integer step, cpt
        integer nx1, nx2, ny1, ny2, nz1, nz2, dnx, dny, dnz
        integer i, j, k, jloc, kloc, npid, nop, nopx, nopy, ierr, nx, ny, nz
        double precision xmx, xmn, ymx, ymn, zmn, zmx, dx, dy, dz
        double precision rho1, rho2, rho3, xmu1, xmu2, xmu3, xmus
        double precision re, re_s, g, st, stv
        double precision sxl, sxr, syl, syr
        double precision sd, stest, stest0, sn, sn0
        double precision adf, adfx, adfy, w
        double precision h_preserve
        common/grid/xmx, xmn, ymx, ymn, zmn, zmx, dx, dy, dz, nx, ny, nz
        common/grid_index/nx1, nx2, ny1, ny2, nz1, nz2, dnx, dny, dnz
        common/mpi/npid, nop, nopx, nopy
        dimension :: sd(nx,ny,nz), stest(nx,ny,nz), stest0(nx,ny,nz)
        dimension :: sn(nx,ny,nz), sn0(nx,ny,nz)
        dimension :: w(nx,ny,nz)
        dimension :: adf(nx,ny,nz), adfx(nx,ny,nz), adfy(nx,ny,nz)
        dimension :: cpt(2)

        adf = 0.d0
        adfx = 0.d0
        adfy = 0.d0
        w = 0.d0
        stest0 = 0.d0
        sn = 0.d0
        sn0 = 0.d0
        ! define thickness of the region outside the interface that is not extrapolated(eg. diffused region)
        ! h_preserve = 0 if only preserve level set inside the interface(sd=0)
        h_preserve = 0.d0
        !--------normal components---------------------------------
        do k = nz1, nz2 
          do j = 2, ny-1
            do i = 2, nx-1

              sxl = sd(i-1,j,1)
              sxr = sd(i+1,j,1)
              syl = sd(i,j-1,1)
              syr = sd(i,j+1,1)
              
              adf(i,j,k) = sqrt( ((sxr-sxl)/2./dx)**2 + ((syr-syl)/2./dy)**2 )

              adfx(i,j,k) = (sxr-sxl)/2./dx / adf(i,j,k)
              adfy(i,j,k) = (syr-syl)/2./dy / adf(i,j,k)

            end do
         end do
       end do
         !--------normal components---------------------------------

         call exchange_u(adfx, adfy, w, cpt)

        !------define sn(diectional derivative of level set s in normal direction)------
           do j = 2, ny-1 
              do i = 2, nx-1
              if( sd(i,j,1).lt.(h_preserve-dx) ) then
              sn(i,j,1) = adfx(i,j,1)*1.d0/2.d0/dx*(stest(i+1,j,1)-stest(i-1,j,1)) + & 
                          adfy(i,j,1)*1.d0/2.d0/dy*(stest(i,j+1,1)-stest(i,j-1,1))
              end if
              end do
           end do
        !------define sn(diectional derivative of level set s in normal direction)------ 


        !------constant extrapolation for sn from inside to outside the interface------ 
        sn0 = sn
        do step = 1, 200

        call levelset_constantprojection(sn0,sd,adfx, adfy, w)
        call exchange_s(sn0, cpt)

        ! retain level set inside solid. Update level set outside solid 
         do j = 2, ny-1
           do i = 2, nx-1
        if(sd(i,j,1).lt.(h_preserve-dx)) then 
        sn0(i,j,1) = sn(i,j,1)
        end if
            end do
         end do
        ! retain level set inside solid. Update level set outside solid

        !call exchange_s(sn0, cpt)
        end do
        sn = sn0
        !------constant extrapolation for sn from inside to outside the interface------ 




        !------linear extrapolation for stest------ 
        !stest0 to store old level set before advection at each time step 
        stest0 = stest
        do step = 1, 200

        call levelset_linearprojection(sn,stest0,sd,adfx, adfy, w)
        call exchange_s(stest0, cpt)

        ! retain level set inside solid. Update level set outside solid 
         do j = ny1-2, ny2+2
           do i = nx1-2, nx2+2
        if(sd(i,j,1).lt.h_preserve) then 
        stest0(i,j,1) = stest(i,j,1)
        end if
            end do
         end do
        ! retain level set inside solid. Update level set outside solid

        end do
        stest = stest0
        !------linear extrapolation for stest------ 


      end subroutine ib_dynamic_grid_advection
