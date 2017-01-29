!!****if* source/Simulation/SimulationMain/unitTest/PoissonWithFFT/Driver_evolveFlash
!!
!! NAME
!!
!!  Driver_evolveFlash
!!
!! SYNOPSIS
!!
!!  Driver_evolveFlash()
!!
!! DESCRIPTION
!!
!!   Solve the possion equation \Delta^2 U(x,y) = F(x,y) on a rectangular grid.
!!   In our case, we want U(x,y) = cos(x+y) to be the solution 
!!   so we take F(x,y) = -2cos(x+y)
!!   
!!   Algorithm:
!!      Define a mesh size, and calculate F(x,y) at mesh points 
!!      Calculate the Fourier coefficients using the cosine basis
!!      Divide the (i,j)'th coefficient by -(i^2+j^2)
!!      Calculate inverse fourier transform to get U(x,y)
!!
!! NOTES
!!     meshsize is decided by the parameters xGridSize and yGridSize
!!     which are integers read from flash.par
!!
!!***

subroutine Driver_evolveFlash()

  use Driver_data, ONLY: dr_globalNumProcs
  use RuntimeParameters_interface, ONLY:  RuntimeParameters_get

  implicit none

#include "constants.h"
#include "Pfftf90.h"   !! Contains PFFT related constants

  !! The size of the grid
  integer :: xGridSize, yGridSize   
  real    :: hx, hy
  real    :: lx,ux,ly,uy
  !! v = values, f = fourier coefficients, 
  !! RHS = Right Hand Side of poisson eqn, Sol = Solution to equation
!!$  real, dimension(:,:) :: vRHS, fRHS, fSol, vSol


  !! pfft related variables
  integer            :: desc 
  integer, parameter :: DIM = 2
  integer, dimension(0:DIM-1) :: size, gsize, psize
  !! size = gridsize, gsize = size + guardcell?, psize = processors on each dimension


  !! local variables
  integer :: i,j
  real    :: inv
  character(len=25) :: fmtstr
  real,dimension(0:31,0:31)::vRHS,fRHS,fSol,vSol

!!$  print *,"Allocating space"
!!$  allocate(vRHS(0:size(0),0:size(1)))
!!$  allocate(fRHS(0:size(0),0:size(1)))
!!$  allocate(fSol(0:size(0),0:size(1)))
!!$  allocate(vSol(0:size(0),0:size(1)))
!!$  print *,"Finished Allocating space"

  !Get the values  of RunTimeParameters
  call RuntimeParameters_get("xGridSize",size(0))
  call RuntimeParameters_get("yGridSize",size(1))
  call RuntimeParameters_get("lx",lx)
  call RuntimeParameters_get("ly",ly)
  call RuntimeParameters_get("ux",ux)
  call RuntimeParameters_get("uy",uy)
  write(*,'("Gridsizes: ",2I4)') size(0), size(1)
  write(*,'("Rectangle: [", 4(F5.2, A))') lx, "-", ux, "] x [", ly, "-", uy, "]"
 
!!$  print *,"Allocating space"
!!$  allocate(vRHS(0:size(0),0:size(1)))
!!$  allocate(fRHS(0:size(0),0:size(1)))
!!$  allocate(fSol(0:size(0),0:size(1)))
!!$  allocate(vSol(0:size(0),0:size(1)))
!!$  print *,"Finished Allocating space"

  print *,"Computing Parameters"
  hx = (ux-lx)/size(0)
  hy = (uy-ly)/size(1)
  
  !! Initialise the RHS matrix
  do i=0, size(0)
     do j=0, size(1)
        vRHS(i,j) = RightHandSide(lx+i*hx,ly+j*hy)
     end do
  end do

  !! compute format string
  write(fmtstr,'("(",I0,"F10.6")')size(0)+1
  print *,"Printing vRHS"
  print fmtstr, vRHS

  !! Parameters for PFFT
  gsize = size
  gsize(0) = gsize(0) + 2 !! No clue why the extra 2? GuardCell?
  !! use one processor along one of the dimensions 
  !! so processor gets access to a whole column of the matrix
  psize(0) = 1
  psize(1) = dr_globalNumProcs
  print *,"Finished Computing Parameters"

  print *,"Calling init_pfft"
  CALL init_pfft(dim=DIM,type=PFFT_COSINE,size=size,gsize=gsize,psize=psize,restore=PFFT_YES,desc=desc)
  print *,"Calling full_transform_forward"
  CALL full_transform_forward(desc,vRHS,fRHS)

  !! print RHS Fourier coefficients
  print *,"Printing fRHS"
  print fmtstr, fRHS

  print *,"Computing F Coeff of Sol"
  !! Compute fourier numbers for Sol
  do i=0, size(0)
     do j=0, size(1)
        inv = i**2 + j**2
        if (inv /= 0) then
           inv = -1/inv
        endif
        fSol(i,j) = fRHS(i,j)*inv
     end do
  end do

  !! print RHS Fourier coefficients
  print *,"Printing fSol"
  print fmtstr, fSol

  print *,"Calling full_transform_inverse"
  CALL full_transform_inverse(desc,fSol,vSol)
  print *,"Printing vSol"
  print fmtstr,vSol

  print *,"Finished calculating Solution"
  print *,"Calling end_pfft"
  CALL end_pfft(desc)
  print *,"Returning from end_pfft"

  return

  contains

  real FUNCTION RightHandSide(x,y)
  REAL, INTENT(IN):: x,y

  RightHandSide = COS(x+y)
  end FUNCTION RightHandSide

end subroutine Driver_evolveFlash


