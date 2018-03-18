subroutine Plasma_velSource(u,v,nrmx,nrmy,s,sigp,ix1,ix2,jy1,jy2,dx,dy)

     use Driver_data, only: dr_nstep, dr_dt

     implicit none

     !Arguements
     real, intent(inout), dimension(:,:,:) :: u,v
     real, intent(inout), dimension(:,:,:) :: sigp
     real, intent(in), dimension(:,:,:) :: s,nrmx,nrmy
     integer, intent(in) :: ix1,ix2,jy1,jy2
     real, intent(in) :: dx,dy

     !Local Variables
     real    :: xS,yS,phiS,phiF
     integer :: k = 1
     integer :: i,j
     real    :: vr = 1.0
     real    :: r,pf,pfx,pfy
     real    :: eps = 1e-13
     real    :: a1x,a2x,a1y,a2y

     phiS = 0.5 - sqrt((xS-0.0)**2+(yS-0.0)**2)

     do j=jy1-1,jy2+1
        do i=ix1-1,ix2+1
  
           pfx = (s(i,j,k)+s(i-1,j,k))/2.0d0
           pfy = (s(i,j,k)+s(i,j-1,k))/2.0d0

           u(i,j,k) = u(i,j,k) + dr_dt*0.5*(nrmx(i,j,k)+nrmx(i-1,j,k))*(vr/(phiS-pfx+eps))
           v(i,j,k) = v(i,j,k) + dr_dt*0.5*(nrmy(i,j,k)+nrmy(i,j-1,k))*(vr/(phiS-pfy+eps))

        end do
     end do


end subroutine Plasma_velSource
