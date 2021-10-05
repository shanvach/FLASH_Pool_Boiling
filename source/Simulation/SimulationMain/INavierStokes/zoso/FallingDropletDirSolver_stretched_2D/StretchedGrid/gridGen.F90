!! NAME
!! -----------------------
!! gridGen
!!
!! DESCRIPTION
!! -----------------------
!! This program generates uniform/stretched grids and
!! save the grid data to an hdf5 file.
!!
!! ARGUMENTS
!!
program gridGen

  implicit none

  integer :: iprocs,jprocs,kprocs
  integer :: nxb,nyb,nzb

  real :: xmin,xmax,ymin,ymax,zmin,zmax
  integer :: ndim,lnx,lny,lnz,nfilter
  real, allocatable, dimension(:) :: strx,stry,strz

  ! Dimensions of the problem
  ndim = 2

  ! Number of elements per processor 
  nxb = 80
  nyb = 80
  nzb = 1
  
  ! Number of processor in each direction
  ! These should be the same as flash.par
  iprocs = 2
  jprocs = 4
  kprocs = 1

  ! Total number of elements in each direction
  lnx = nxb*iprocs
  lny = nyb*jprocs
  lnz = nzb*kprocs

  ! Physical domain lengths
  ! These should be the same as flash.par
  xmin = -2.0
  xmax =  2.0
  ymin =  0.0
  ymax =  8.0
  zmin = -1.0
  zmax =  1.0

  ! Number of filtering operation
  nfilter = 20

  allocate(strx(iprocs),stry(jprocs),strz(kprocs))

  ! Stretching coefficients.
  ! The stretching values better be close to 1.
  ! Otherwise conservation errors may increase.
  ! If the stretching values are set much different than 1 then
  ! number of filtering operations can be increased to stabilize it.

  ! In order to decrease the cell spacing towards the end of the block,
  ! stretching coefficient should be less then 1, and bigger than 1 vice versa.

  ! Number of stretching coefficients should be equal to number of
  ! processors in respected direction.
  strx = (/0.99,1.01/)
  stry = (/1.0,1.0,1.0,1.0/)
  strz = (/1.0,0.5,1.5,1.0/)

  if(ndim==3 .and. nzb==1) stop "It is a 3D problem. Check configuration in z."

  call gridGen_main(nxb,nyb,nzb,lnx,lny,lnz,ndim,iprocs,jprocs,kprocs,&
                    xmin,xmax,ymin,ymax,zmin,zmax,strx,stry,strz,nfilter)
 
  deallocate(strx,stry,strz)

end program

subroutine gridGen_main(nxb,nyb,nzb,lnx,lny,lnz,ndim,iprocs,jprocs,kprocs,&
                        xmin,xmax,ymin,ymax,zmin,zmax,strx,stry,strz,nfilter)
 
  use HDF5

  implicit none

  integer, intent(in) :: nxb,nyb,nzb,lnx,lny,lnz,ndim
  integer, intent(in) :: iprocs,jprocs,kprocs,nfilter
  real, intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax
 
  character(len=11), parameter :: filename = "initGrid.h5"
  character(len=7) :: dset_name
  integer :: error             ! Error code
  integer :: dset_rank   ! Rank of the data set

  integer(HSIZE_T), dimension(2) :: dset_dims   ! Dimensions of the data set
  integer(HID_T) :: dset_id     ! Dataset identifier
  integer(HID_T) :: file_id     ! File identifier
  integer(HID_T) :: dset_space  ! Dataspace identifier

  double precision, dimension(2,nxb) :: x_coordsLoc
  double precision, dimension(lnx+1) :: x_temp
  double precision, dimension(lnx+1) :: x_Faces

  double precision, dimension(2,nyb) :: y_coordsLoc
  double precision, dimension(lny+1) :: y_temp
  double precision, dimension(lny+1) :: y_Faces

  double precision, dimension(2,nzb) :: z_coordsLoc
  double precision, dimension(lnz+1) :: z_temp
  double precision, dimension(lnz+1) :: z_Faces

  real, dimension(iprocs), intent(in) :: strx
  real, dimension(jprocs), intent(in) :: stry 
  real, dimension(jprocs), intent(in) :: strz 
 
  real :: xmin_loc,xmax_loc,str
  real :: ymin_loc,ymax_loc
  real :: zmin_loc,zmax_loc

  integer :: i,j,k

  x_coordsLoc = 0.0
  x_temp = 0.0

!!  =================================================
!!  x-Face Coordinates 
!!  =================================================

  do i=1,iprocs
     str = strx(i)
     xmin_loc = (i-1)*((xmax-xmin)/iprocs)+xmin
     xmax_loc = (i  )*((xmax-xmin)/iprocs)+xmin     
     call createDomain_geomSeries(xmin_loc,xmax_loc,str,nxb,x_coordsLoc)
     if(str .lt. 1.0) then
       x_Faces((i-1)*nxb+2:i*nxb) = x_coordsLoc(2,:)
       x_Faces((i-1)*nxb+1) = xmin_loc
     else
       x_Faces((i-1)*nxb+1:i*nxb) = x_coordsLoc(1,:)
       x_Faces((i)*nxb+1) = xmax_loc
     endif 
  enddo

  !Filter/Smooth grid
  x_temp = x_Faces
  do j=1,nfilter
    do i=2,lnx
      x_temp(i)=0.25*(x_Faces(i-1)+2.*x_Faces(i)+x_Faces(i+1))
    enddo
    x_Faces = x_temp
  enddo

  dset_name = "/xFaces"
  dset_rank = 2
  dset_dims = (/lnx+1,1/)
 
  call h5open_f(error)

  ! Create a new file
  call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,error)

  ! Create the data space for the  dataset. 
  call h5screate_simple_f(dset_rank,dset_dims,dset_space,error)

  ! Create the dataset
  call h5dcreate_f(file_id,dset_name,H5T_NATIVE_DOUBLE,dset_space,dset_id,error)

  call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,x_Faces(:),dset_dims,error)

!!  =================================================
!!  y-Face Coordinates 
!!  =================================================

  if(ndim == 2) then

  do j=1,jprocs
     str = stry(j)
     ymin_loc = (j-1)*((ymax-ymin)/jprocs)+ymin
     ymax_loc = (j  )*((ymax-ymin)/jprocs)+ymin     
     call createDomain_geomSeries(ymin_loc,ymax_loc,str,nyb,y_coordsLoc)
     if(str .lt. 1.0) then
       !Populate the first location
       y_Faces((j-1)*nyb+2:j*nyb) = y_coordsLoc(2,:)
       !Populate the last location
       y_Faces((j-1)*nyb+1) = ymin_loc
     else
       y_Faces((j-1)*nyb+1:j*nyb) = y_coordsLoc(1,:)
       y_Faces((j)*nyb+1) = ymax_loc
     endif
  enddo

  !Filter/Smooth grid
  y_temp = y_Faces
  do j=1,nfilter
    do i=2,lny
      y_temp(i)=0.25*(y_Faces(i-1)+2.*y_Faces(i)+y_Faces(i+1))
    enddo
    y_Faces = y_temp
  enddo

  dset_name = "/yFaces"
  dset_rank = 2
  dset_dims = (/lny+1,1/)
   
  ! Create the data space for the  dataset. 
  call h5screate_simple_f(dset_rank,dset_dims,dset_space,error)

  ! Create the dataset
  call h5dcreate_f(file_id,dset_name,H5T_NATIVE_DOUBLE,dset_space,dset_id,error)

  call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,y_Faces(:),dset_dims,error)

  endif

!!  =================================================
!!  z-Face Coordinates 
!!  =================================================

  if(ndim==3) then

  do k=1,jprocs
     str = strz(k)
     zmin_loc = (k-1)*((zmax-zmin)/kprocs)+zmin
     zmax_loc = (k  )*((zmax-zmin)/kprocs)+zmin     
     call createDomain_geomSeries(zmin_loc,zmax_loc,str,nzb,z_coordsLoc)
     if(str .lt. 1.0) then
       z_Faces((k-1)*nzb+2:k*nzb) = z_coordsLoc(2,:)
       z_Faces((k-1)*nzb+1) = zmin_loc
     else
       z_Faces((k-1)*nzb+1:k*nzb) = z_coordsLoc(1,:)
       z_Faces((k)*nzb+1) = zmax_loc
     endif
  enddo

  !Filter/Smooth grid
  z_temp = z_Faces
  do j=1,nfilter
    do i=2,lnz
      z_temp(i)=0.25*(z_Faces(i-1)+2.*z_Faces(i)+z_Faces(i+1))
    enddo
    z_Faces = z_temp
  enddo

  dset_name = "/zFaces"
  dset_rank = 2
  dset_dims = (/lnz+1,1/)
   
  ! Create the data space for the  dataset. 
  call h5screate_simple_f(dset_rank,dset_dims,dset_space,error)

  ! Create the dataset
  call h5dcreate_f(file_id,dset_name,H5T_NATIVE_DOUBLE,dset_space,dset_id,error)

  call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,z_Faces(:),dset_dims,error)

  endif

  ! Close the dataspace for the dataset.
  call h5sclose_f(dset_space, error)

  ! Close the dataset.
  call h5dclose_f(dset_id, error)

  ! Close the file.
  call h5fclose_f(file_id, error)

  ! Close FORTRAN interface.
  call h5close_f(error) 

end subroutine

subroutine createDomain_geomSeries(xmin,xmax,str,lnx,x_coords)
  
  implicit none
  
  real, intent(in) :: xmin,xmax,str
  integer, intent(in) :: lnx
  double precision, dimension(2,lnx), intent(inout) :: x_coords

  real :: s,dx
  integer :: i

  s = 0.0
  x_coords = 0.0
 
  do i=1,lnx
     s = s + str**(i-1)
  enddo

  dx = (xmax-xmin)/s

  x_coords(1,1) = xmin
  x_coords(2,1) = xmin+dx
 
  do i=2,lnx
     x_coords(1,i) = x_coords(1,i-1) + dx*(str**(i-1))
     x_coords(2,i) = x_coords(2,i-1) + dx*(str**(i-1))
  enddo

end subroutine
