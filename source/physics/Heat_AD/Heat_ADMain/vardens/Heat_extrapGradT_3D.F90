subroutine Heat_extrapGradT_3D(Tnl,Tnv,T,s,pf,dx,dy,dz,nx,ny,nz,ix1,ix2,jy1,jy2,kz1,kz2,Tnl_res,Tnv_res,mflg)

#include "Flash.h"

    implicit none
    real, dimension(:,:,:), intent(inout) :: Tnl,Tnv
    real, dimension(:,:,:), intent(in) :: T,s,pf,nx,ny,nz,mflg
    real, intent(in) :: dx,dy,dz
    integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
    real, intent(out) :: Tnl_res,Tnv_res

    integer :: i,j,k

    real, dimension(NXB,NYB,NZB) :: nxconv,nyconv,nzconv,nx_plus,nx_mins,ny_plus,ny_mins,nz_plus,nz_mins
    real, dimension(NXB,NYB,NZB) :: Tlx_plus, Tlx_mins, Tly_plus, Tly_mins, Tlz_plus, Tlz_mins
    real, dimension(NXB,NYB,NZB) :: Tvx_plus, Tvx_mins, Tvy_plus, Tvy_mins, Tvz_plus, Tvz_mins

    real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: Tnl_o,Tnv_o,Tnl_i,Tnv_i

    real :: dt_ext,maxN

    integer :: l_counter, v_counter, step

    l_counter = 0
    v_counter = 0

    maxN = max(maxval(abs(nx)),maxval(abs(ny)))

    dt_ext = 0.5d0*(dx/maxN)
    !dt_ext = 1E-6
    !dt_ext = 0.5d0/(1.0d0/dx + 1.0d0/dy)
    !dt_ext = dx

    Tnl_i = Tnl
    Tnv_i = Tnv

    Tnl_o = Tnl

    Tnv_o = Tnv

    nxconv = sign(1.,s(ix1:ix2,jy1:jy2,kz1:kz2))*nx(ix1:ix2,jy1:jy2,kz1:kz2)
    nyconv = sign(1.,s(ix1:ix2,jy1:jy2,kz1:kz2))*ny(ix1:ix2,jy1:jy2,kz1:kz2)
    nzconv = sign(1.,s(ix1:ix2,jy1:jy2,kz1:kz2))*nz(ix1:ix2,jy1:jy2,kz1:kz2)

    nx_plus = max(nxconv,0.)
    nx_mins = min(nxconv,0.)

    ny_plus = max(nyconv,0.)
    ny_mins = min(nyconv,0.)

    nz_plus = max(nzconv,0.)
    nz_mins = min(nzconv,0.)

    Tlx_plus = (Tnl_o(ix1+1:ix2+1,jy1:jy2,kz1:kz2)-Tnl_o(ix1:ix2,jy1:jy2,kz1:kz2))/dx
    Tlx_mins = (Tnl_o(ix1:ix2,jy1:jy2,kz1:kz2)-Tnl_o(ix1-1:ix2-1,jy1:jy2,kz1:kz2))/dx

    Tly_plus = (Tnl_o(ix1:ix2,jy1+1:jy2+1,kz1:kz2)-Tnl_o(ix1:ix2,jy1:jy2,kz1:kz2))/dy
    Tly_mins = (Tnl_o(ix1:ix2,jy1:jy2,kz1:kz2)-Tnl_o(ix1:ix2,jy1-1:jy2-1,kz1:kz2))/dy

    Tlz_plus = (Tnl_o(ix1:ix2,jy1:jy2,kz1+1:kz2+1)-Tnl_o(ix1:ix2,jy1:jy2,kz1:kz2))/dz
    Tlz_mins = (Tnl_o(ix1:ix2,jy1:jy2,kz1:kz2)-Tnl_o(ix1:ix2,jy1:jy2,kz1-1:kz2-1))/dz

    Tvx_plus = (Tnv_o(ix1+1:ix2+1,jy1:jy2,kz1:kz2)-Tnv_o(ix1:ix2,jy1:jy2,kz1:kz2))/dx
    Tvx_mins = (Tnv_o(ix1:ix2,jy1:jy2,kz1:kz2)-Tnv_o(ix1-1:ix2-1,jy1:jy2,kz1:kz2))/dx

    Tvy_plus = (Tnv_o(ix1:ix2,jy1+1:jy2+1,kz1:kz2)-Tnv_o(ix1:ix2,jy1:jy2,kz1:kz2))/dy
    Tvy_mins = (Tnv_o(ix1:ix2,jy1:jy2,kz1:kz2)-Tnv_o(ix1:ix2,jy1-1:jy2-1,kz1:kz2))/dy

    Tvz_plus = (Tnv_o(ix1:ix2,jy1:jy2,kz1+1:kz2+1)-Tnv_o(ix1:ix2,jy1:jy2,kz1:kz2))/dz
    Tvz_mins = (Tnv_o(ix1:ix2,jy1:jy2,kz1:kz2)-Tnv_o(ix1:ix2,jy1:jy2,kz1-1:kz2-1))/dz

    !Tnl(ix1:ix2,jy1:jy2,kz1:kz2) = Tnl_i(ix1:ix2,jy1:jy2,kz1:kz2) + dt_ext*pf(ix1:ix2,jy1:jy2,kz1:kz2)*&
    !                                                              (-nx_mins*Tlx_plus-nx_plus*Tlx_mins &
    !                                                               -ny_mins*Tly_plus-ny_plus*Tly_mins &
    !                                                               -nz_mins*Tlz_plus-nz_plus*Tlz_mins)


    !Tnv(ix1:ix2,jy1:jy2,kz1:kz2) = Tnv_i(ix1:ix2,jy1:jy2,kz1:kz2) + dt_ext*(1.0-pf(ix1:ix2,jy1:jy2,kz1:kz2))*&
    !                                                              (-nx_mins*Tvx_plus-nx_plus*Tvx_mins &
    !                                                               -ny_mins*Tvy_plus-ny_plus*Tvy_mins &
    !                                                               -nz_mins*Tvz_plus-nz_plus*Tvz_mins)

  do k =kz1,kz2
   do j = jy1,jy2
     do i = ix1,ix2

       !if((s(i,j,k)*s(i+1,j,k) .le. 0.) .or. &
       !   (s(i,j,k)*s(i-1,j,k) .le. 0.) .or. &
       !   (s(i,j,k)*s(i,j+1,k) .le. 0.) .or. &
       !   (s(i,j,k)*s(i,j-1,k) .le. 0.) .or. &
       !   (s(i,j,k)*s(i,j,k+1) .le. 0.) .or. &
       !   (s(i,j,k)*s(i,j,k-1) .le. 0.)) then

        Tnl(i,j,k) = Tnl_i(i,j,k) + mflg(i,j,k)*dt_ext*pf(i,j,k)*(-nx_mins(i-ix1+1,j-jy1+1,k-kz1+1)*Tlx_plus(i-ix1+1,j-jy1+1,k-kz1+1) &
                                                                  -nx_plus(i-ix1+1,j-jy1+1,k-kz1+1)*Tlx_mins(i-ix1+1,j-jy1+1,k-kz1+1) &
                                                                  -ny_mins(i-ix1+1,j-jy1+1,k-kz1+1)*Tly_plus(i-ix1+1,j-jy1+1,k-kz1+1) &
                                                                  -ny_plus(i-ix1+1,j-jy1+1,k-kz1+1)*Tly_mins(i-ix1+1,j-jy1+1,k-kz1+1) &
                                                                  -nz_mins(i-ix1+1,j-jy1+1,k-kz1+1)*Tlz_plus(i-ix1+1,j-jy1+1,k-kz1+1) &
                                                                  -nz_plus(i-ix1+1,j-jy1+1,k-kz1+1)*Tlz_mins(i-ix1+1,j-jy1+1,k-kz1+1))

        Tnv(i,j,k) = Tnv_i(i,j,k) + mflg(i,j,k)*dt_ext*(1.0-pf(i,j,k))*(-nx_mins(i-ix1+1,j-jy1+1,k-kz1+1)*Tvx_plus(i-ix1+1,j-jy1+1,k-kz1+1) &
                                                                        -nx_plus(i-ix1+1,j-jy1+1,k-kz1+1)*Tvx_mins(i-ix1+1,j-jy1+1,k-kz1+1) &
                                                                        -ny_mins(i-ix1+1,j-jy1+1,k-kz1+1)*Tvy_plus(i-ix1+1,j-jy1+1,k-kz1+1) &
                                                                        -ny_plus(i-ix1+1,j-jy1+1,k-kz1+1)*Tvy_mins(i-ix1+1,j-jy1+1,k-kz1+1) &
                                                                        -nz_mins(i-ix1+1,j-jy1+1,k-kz1+1)*Tvz_plus(i-ix1+1,j-jy1+1,k-kz1+1) &
                                                                        -nz_plus(i-ix1+1,j-jy1+1,k-kz1+1)*Tvz_mins(i-ix1+1,j-jy1+1,k-kz1+1))

       !endif

     end do
  end do
 end do

    !do k =kz1,kz2
    ! do j =jy1,jy2
    !  do i=ix1,ix2

    !  if((s(i,j,k)*s(i+1,j,k) .le. 0.) .or. &
    !     (s(i,j,k)*s(i-1,j,k) .le. 0.) .or. &
    !     (s(i,j,k)*s(i,j+1,k) .le. 0.) .or. &
    !     (s(i,j,k)*s(i,j-1,k) .le. 0.) .or. &
    !     (s(i,j,k)*s(i,j,k+1) .le. 0.) .or. &
    !     (s(i,j,k)*s(i,j,k-1) .le. 0.)) then

    !     Tnl_res = Tnl_res + (Tnl_o(i,j,k)-Tnl(i,j,k))**2
    !     Tnv_res = Tnv_res + (Tnv_o(i,j,k)-Tnv(i,j,k))**2

    !  end if


    !  end do
    ! end do
    !end do

    Tnl_res = sum(sum(sum((Tnl_o(:,:,:)-Tnl(:,:,:))**2,1),1))
    Tnv_res = sum(sum(sum((Tnv_o(:,:,:)-Tnv(:,:,:))**2,1),1))

    Tnl_res = sqrt(Tnl_res/size(Tnl))
    Tnv_res = sqrt(Tnv_res/size(Tnv))

end subroutine Heat_extrapGradT_3D
