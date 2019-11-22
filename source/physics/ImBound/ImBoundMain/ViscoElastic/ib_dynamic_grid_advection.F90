!===============================================================================
!!
!! subroutine ib_dynamic_grid_advection 
!!
!! extrapolating a scalar field from inside to the outside of a region defined by sd
!!
!! THIS NEEDS TO BE BROKEN DOWN INTO MULTIPLE SUBROUTINES IN ACCORDANCE
!  WITH imBound.F90
!===============================================================================



      !subroutine ib_dynamic_grid_advection(sn, sd, stest, adfx, adfy, cpt)
      subroutine ib_dynamic_grid_advection(sn,sd,stest,adf,adfx,adfy,&
                                           cpt,dx,dy,dz,ix1,ix2,jy1,jy2)
        implicit none
        include 'mpif.h' !where is this?

        real, dimension(:,:,:), intent(inout) :: sn,adf,adfx,adfy
        real, dimension(:,:,:), intent(in)    :: sd,stest
        integer, dimension(:),  intent(in)    :: cpt

        real, intent(in)    :: dx,dy,dz
        integer, intent(in) :: ix1,ix2,jy1,jy2
        !end interface header

        integer :: step,i,j,k

        real :: sxl, sxr, syl, syr
        !end buit-in header

        !functions called
        !exchange_u(adfx, adfy, w, cpt)
        !call levelset_constantprojection(sn0,sd,adfx, adfy, w)
        !call exchange_s(sn0, cpt)
        !call levelset_linearprojection(sn,stest0,sd,adfx, adfy, w)
        !call exchange_s(stest0, cpt)
        
        !ib_dynamic_grid_advection(**args)
        !double precision sn, sd, stest, adfx, adfy
        !integer cpt
        !vars that are common/grid
        !common/grid/xmx, xmn, ymx, ymn, zmn, zmx, dx, dy, dz, nx, ny, nz
        !double precision xmx, xmn, ymx, ymn, zmn, zmx, dx, dy, dz, nx, ny, nz
        !vars that are common grid_index
        !common/grid_index/nx1, nx2, ny1, ny2, nz1, nz2, dnx, dny, dnz
        !integer :: nx1, nx2, ny1, ny2, nz1, nz2, dnx, dny, dnz
        !vars that are common to common/mpi
        !common/mpi/npid, nop, nopx, nopy
        !integer :: npid, nop, nopx, nopy
        !vars created in this function
        !integer :: i, j, k, jloc, kloc, ierr
        !double precision rho1, rho2, rho3, xmu1, xmu2, xmu3, xmus
        !double precision re, re_s, g, st, stv
        !double precision sxl, sxr, syl, syr
        !double precision stest0, sn0
        !double precision adf, w
        !double precision h_preserve
        !dimension :: stest0(nx,ny,nz), sn0(nx,ny,nz)
        !dimension :: adf(nx,ny,nz), w(nx,ny,nz)

        !integer :: nx1, nx2, ny1, ny2, nz1, nz2, dnx, dny, dnz
        !integer :: i, j, k, jloc, kloc, ierr

        !double precision xmx, xmn, ymx, ymn, zmn, zmx, dx, dy, dz

        !double precision rho1, rho2, rho3, xmu1, xmu2, xmu3, xmus
        !double precision re, re_s, g, st, stv
        !double precision sxl, sxr, syl, syr
        !double precision stest, stest0, sn0
        !double precision adf, w
        !double precision h_preserve

        !these need to be included from a data file not directly here
        !common/grid/xmx, xmn, ymx, ymn, zmn, zmx, dx, dy, dz, nx, ny, nz
        !common/grid_index/nx1, nx2, ny1, ny2, nz1, nz2, dnx, dny, dnz
        !common/mpi/npid, nop, nopx, nopy

        !dimension :: stest0(nx,ny,nz), sn0(nx,ny,nz)
        !dimension :: w(nx,ny,nz)
        !dimension :: adf(nx,ny,nz)

        !integer step!, cpt
        !dimension :: sd(nx,ny,nz), stest(nx,ny,nz),
        !dimension :: sn(nx,ny,nz), adfx(nx,ny,nz), adfy(nx,ny,nz)
        !dimension :: cpt(2)

        !**args
        adfx = 0.d0
        adfy = 0.d0
        sn = 0.d0

        !created in this function, now passed as args
        adf = 0.d0
        w = 0.d0
        stest0 = 0.d0
        sn0 = 0.d0
        ! define thickness of the region outside the interface that is not extrapolated(eg. diffused region)
        ! h_preserve = 0 if only preserve level set inside the interface(sd=0)
        h_preserve = 0.d0

        !this obtains normal components (nx,ny) of interface
        !FLASH has a method to do this, so we may not need this next
        !loop
        !--------normal components---------------------------------
        !do k = nz1,nz2 
        k = 1
          !do j = 2, ny-1
          do j = jy1,jy2
            !do i = 2, nx-1
            do i = ix1,ix2

              sxl = sd(i-1,j,1)
              sxr = sd(i+1,j,1)
              syl = sd(i,j-1,1)
              syr = sd(i,j+1,1)
              
              adf(i,j,k) = sqrt( ((sxr-sxl)/2./dx)**2 + ((syr-syl)/2./dy)**2 )

              adfx(i,j,k) = (sxr-sxl)/2./dx / adf(i,j,k)
              adfy(i,j,k) = (syr-syl)/2./dy / adf(i,j,k)

            end do
         end do
       !end do
         !--------normal components---------------------------------

         call exchange_u(adfx, adfy, w, cpt)
         !This is accomplished by filling guard cells, so need to break
         !up this into multiple subroutines in order to fill guard cells
    
       
        !------define sn(diectional derivative of level set s in normal direction)------
           !do j = 2, ny-1 
           do j = jy1,jy2
              !do i = 2, nx-1
              do i = ix1,ix2
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
        !This is accomplished by filling guard cells, so need to break
         !up this into multiple subroutines in order to fill guard cells

        ! retain level set inside solid. Update level set outside solid 
         !do j = 2, ny-1
         do j = jy1,jy2
           !do i = 2, nx-1
           do i = ix1,ix2
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
        !This is accomplished by filling guard cells, so need to break
         !up this into multiple subroutines in order to fill guard cells

        ! retain level set inside solid. Update level set outside solid 
         !do j = ny1-2, ny2+2
         do j = jy1,jy2
           !do i = nx1-2, nx2+2
           do i = ix1,ix2
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
