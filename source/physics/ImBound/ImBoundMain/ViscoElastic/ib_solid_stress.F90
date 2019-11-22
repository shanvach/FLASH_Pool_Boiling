!===============================================================================
!!
!! subroutine ib_solid_stress
!!
!! calculate deformation gradient tensor , solid strain tensor and strain magnitude 
!!
!! args inout:: Tau       - strain tensor (2D=4elements)
!! args in   :: sd        - level set for phase interface
!!           :: sX,sY     - dynamic grid (X,Y)
!!           :: A,AT,A_inv- deformation gradient tensor, transpose, and
!!                          inverse
!!           :: Taux,Tauy - Div(Tau),divergence of tensor 
!!           :: Taum      - magnitude(Taux,Tauy)
!!           :: dx,dy     - spacing in x and y directions
!!           :: ix1,ix2   - low,high block x indeces
!!           :: jy1,jy2   - low,high block y indeces 
!===============================================================================
        subroutine ib_solid_stress(sd,sY,sX,A,AT,A_inv,&
                                   Taux,Tauy,Taum,Tau,&
                                   dx,dy,ix1,ix2,jy1,jy2)   
        implicit none
        !include 'mpif.h'
        real, dimension(:,:,:), intent(inout) :: Tau,Taux,Tauy,Taum
        real, dimension(:,:,:), intent(in)    :: sd,sY,sX,A,AT,A_inv

        real, intent(in)    :: dx,dy
        integer, intent(in) :: ix1,ix2,jy1,jy2

        integer :: i,j,k
        real    :: ul,ur,vl,vr
        !end header

        !vars created in this subroutine
        !A_inv(nx,ny,4), A(nx,ny,4), AT(nx,ny,4)

        !integer i, j, k, nx, ny, nz
        !integer nx1, nx2, ny1, ny2, nz1, nz2
        !double precision Taux, Tauy, Tau_magnitude
        !double precision A_inv, A, AT, Tau
        !double precision sY, sX, sd
        !double precision dx, dy, dz
        !dimension   ::   sY(nx,ny,nz), sX(nx,ny,nz), sd(nx,ny,nz)
        !dimension   ::   A_inv(nx,ny,4), A(nx,ny,4), AT(nx,ny,4), Tau(nx,ny,4)
        !dimension   ::   Taux(nx,ny,nz), Tauy(nx,ny,nz), Tau_magnitude(nx,ny,nz)
        !common/grid/     dx, dy, dz, nx, ny, nz
        !common/grid_index/nx1, nx2, ny1, ny2, nz1, nz2

        A_inv = 0.d0
        A     = 0.d0
        AT    = 0.d0

        Tau  = 0.d0
        Taux = 0.d0
        Tauy = 0.d0
        Taum = 0.d0 

        !A(i,j,ind) is the deformation gradient tensor in 2D 
        !ind is the index of each component in the deformation gradient tensor. dX/dx:1; dX/dy:2; dY/dx:3; dY/dy:4.
        k = 1
        do j = jy1,jy2
          do i = ix1,ix2

          A_inv(i,j,1) = 1.d0/(2*dx)*(sX(i+1,j,k)-sX(i-1,j,k))
          A_inv(i,j,2) = 1.d0/(2*dy)*(sX(i,j+1,k)-sX(i,j-1,k))
          A_inv(i,j,3) = 1.d0/(2*dx)*(sY(i+1,j,k)-sY(i-1,j,k))
          A_inv(i,j,4) = 1.d0/(2*dy)*(sY(i,j+1,k)-sY(i,j-1,k))


          if(abs(A_inv(i,j,1)*A_inv(i,j,4)-A_inv(i,j,2)*A_inv(i,j,3)).gt.1E-12) then
           A(i,j,1) = 1.d0/(A_inv(i,j,1)*A_inv(i,j,4)-A_inv(i,j,2)*A_inv(i,j,3))  &
                      *A_inv(i,j,4)
           A(i,j,2) = 1.d0/(A_inv(i,j,1)*A_inv(i,j,4)-A_inv(i,j,2)*A_inv(i,j,3))  &
                      *(-A_inv(i,j,2)) 
           A(i,j,3) = 1.d0/(A_inv(i,j,1)*A_inv(i,j,4)-A_inv(i,j,2)*A_inv(i,j,3))  &
                      *(-A_inv(i,j,3)) 
           A(i,j,4) = 1.d0/(A_inv(i,j,1)*A_inv(i,j,4)-A_inv(i,j,2)*A_inv(i,j,3))  &
                      *A_inv(i,j,1) 
          end if

          AT(i,j,1) = A(i,j,1)
          AT(i,j,2) = A(i,j,3)
          AT(i,j,3) = A(i,j,2)
          AT(i,j,4) = A(i,j,4)

          !Tau is the strain tensor
          Tau(i,j,1) = (A(i,j,1)*AT(i,j,1)+A(i,j,2)*AT(i,j,3)-1)
          Tau(i,j,2) = (A(i,j,1)*AT(i,j,2)+A(i,j,2)*AT(i,j,4))
          Tau(i,j,3) = (A(i,j,3)*AT(i,j,1)+A(i,j,4)*AT(i,j,3))
          Tau(i,j,4) = (A(i,j,3)*AT(i,j,2)+A(i,j,4)*AT(i,j,4)-1)
        
          end do
        end do

        !divergence of Tau gives solid stress in x and y directions
        k = 1
        do j = jy1,jy2
          do i = ix1,ix2
          !if(sd(i,j,1).le.0.d0) then
          Taux(i,j,k) = 1.d0/(2*dx)*(Tau(i+1,j,1)-Tau(i-1,j,1)) + 1.d0/(2*dy)*(Tau(i,j+1,2)-Tau(i,j-1,2)) 
          Tauy(i,j,k) = 1.d0/(2*dx)*(Tau(i+1,j,3)-Tau(i-1,j,3)) + 1.d0/(2*dy)*(Tau(i,j+1,4)-Tau(i,j-1,4))
          Taum(i,j,k) = sqrt(Taux(i,j,k)**2 + Tauy(i,j,k)**2)
        !end if 
           end do
        end do

      end subroutine ib_solid_stress
