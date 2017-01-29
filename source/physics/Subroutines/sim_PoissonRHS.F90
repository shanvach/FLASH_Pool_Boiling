subroutine sim_PoissonRHS(prhs,phat,p,rho1x,rho2x,rho1y,rho2y,rho1z,rho2z,&
                          rho,ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz)

  !____________________ARGUMENTS________________________!

   ! prhs - variable containing value of p from n-1 time level and storing result
   ! phat - variable to store value of phat
   ! rho1x,rho2x - phase densities in x direction
   ! rho1y,rho2y - phase densities in y direction
   ! rho1z,rho2z - phase densities in z direction
   ! rho - minimum of (rho1/rho2) and (rho2/rho2)
   ! ix1,ix2 - block limits in x direction
   ! jy1,jy2 - block limits in y direction
   ! kz1,kz2 - block limits in k direction
   ! dx,dy,dz - block deltas

   implicit none

   real,dimension(:,:,:),intent(inout) :: phat,prhs
   real,dimension(:.:.:),intent(in) :: p,rho1x,rho2x,rho1y,rho2y,rho1z,rho2z
   integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
   real, intent(in) :: dx,dy,dz
   real, intent(in) :: rho


  !____________________LOCAL VARIABLES__________________!

#if NDIM == 2
   real,dimension(NXB+2*NGUARD,NYB+2*NGUARD,1) :: ptemp
#endif

#ifdef NDIM == 3
   real,dimension(NXB+2*NGUARD,NYB+2*NGUARD,NZB+2*NGUARD) :: ptemp
#endif

   real :: pxplus,pxmins,pyplus,pymins,pzplus,pzmins
   real :: Mdxplus,Mdxmins,Mdyplus,Mdymins,Mdzplus,Mdzmins

   integer :: i,j,k


  !_____________________FUNCTION BODY__________________!

   phat = 2*p - prhs


#if NDIM == 2
   k = 1

   do j=jy1,jy2     
      do i=ix1,ix2

         Mdyplus = (rho1y(i,j+1,k) + rho2y(i,j+1,k))
         Mdxplus = (rho1x(i+1,j,k) + rho2x(i+1,j,k))

         Mdymins = (rho1y(i,j,k) + rho2y(i,j,k))
         Mdxmins = (rho1x(i,j,k) + rho2x(i,j,k))

         pxplus = (1 - (rho*Mdxplus))*((phat(i+1,j,k)-phat(i,j,k))/dx)
         pyplus = (1 - (rho*Mdyplus))*((phat(i,j+1,k)-phat(i,j,k))/dy)

         pxmins = (1 - (rho*Mdxmins))*((phat(i,j,k)-phat(i-1,j,k))/dx)
         pymins = (1 - (rho*Mdymins))*((phat(i,j,k)-phat(i,j-1,k))/dy)


         ptemp(i,j,k) = (pxplus-pxmins)/dx + (pyplus-pymins)/dy

      end do
   end do

   prhs = ptemp

#endif

#ifdef NDIM == 3

 do k=kz1,kz2
   do j=jy1,jy2
      do i=ix1,ix2

         Mdyplus = (rho1y(i,j+1,k) + rho2y(i,j+1,k))
         Mdxplus = (rho1x(i+1,j,k) + rho2x(i+1,j,k))
         Mdzplus = (rho1z(i,j,k+1) + rho2z(i,j,k+1))

         Mdymins = (rho1y(i,j,k) + rho2y(i,j,k))
         Mdxmins = (rho1x(i,j,k) + rho2x(i,j,k))
         Mdzmins = (rho1z(i,j,k) + rho2z(i,j,k))

         pxplus = (1 - (rho*Mdxplus))*((phat(i+1,j,k)-phat(i,j,k))/dx)
         pyplus = (1 - (rho*Mdyplus))*((phat(i,j+1,k)-phat(i,j,k))/dy)
         pzplus = (1 - (rho*Mdzplus))*((phat(i,j,k+1)-phat(i,j,k))/dz)

         pxmins = (1 - (rho*Mdxmins))*((phat(i,j,k)-phat(i-1,j,k))/dx)
         pymins = (1 - (rho*Mdymins))*((phat(i,j,k)-phat(i,j-1,k))/dy)
         pzmins = (1 - (rho*Mdzmins))*((phat(i,j,k)-phat(i,j,k-1))/dz)


         ptemp(i,j,k) = (pxplus-pxmins)/dx + (pyplus-pymins)/dy + (pzplus-pzmins)/dz

      end do
   end do
 end do

   prhs = ptemp

#endif

end subroutine sim_PoissonRHS
