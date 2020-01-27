!!****if* source/Grid/GridSolvers/Pfft/RegularGridSolver/gr_pfftPoissonDirect
!!
!! NAME
!!
!!  gr_pfftPoissonDirect
!!
!!
!! SYNOPSIS
!!
!!  call gr_pfftPoissonDirect(integer(IN)                :: iDirection,
!!                            integer(IN)                :: solveflag,
!!                            integer(IN)                :: inSize,
!!                            integer(IN)                :: localSize(MDIM),
!!                            integer(IN)                :: globalSize(MDIM),
!!                            integer(IN)                :: transformType(MDIM),
!!                            real(IN) ,dimension(inSize):: inArray(:),
!!                            real(OUT),dimension(inSize):: outArray(:))
!!
!! DESCRIPTION
!!
!!   Poisson solver routine.  This module implements the direct based method 
!!   for periodic, neumann, and dirichlet problems on a regular grid. 
!!   Mixed boundaries for a given axis are not supported.
!!
!!
!! ARGUMENTS
!!
!!   iDirection  - direction of the transform, valid values are 
!!                 PFFT_FORWARD and PFFT_INVERSE 
!!   solveflag   - Indicates the solvers to be applied.
!!                 solveflag==TRIG_DRCT      => Transform in X, Direct Solve in Y
!!                 solveflag==TRIG_TRIG_DRCT => Transform in X, Direct Solve in Y, Z
!!                 solveflag==TRIG_DRCT_DRCT => Transform in X and Y, Direct Solve in Z
!!   inSize      - size of inArray and outArray
!!   localSize   - the local bounds (on myPe) of the original 3D data to be transformed
!!   globalSize  - global size of the 3D data to be transformed
!!   transformType - The type if transform to be applied along each of the dimensions
!!   inArray       - single dimension array containing linearized data to be transformed
!!   outArray      - single dimension array containing transformed linearized data
!!  
!!***

subroutine gr_pfftPoissonDirect (iDirection, solveflag, inSize, localSize, globalSize, transformType, inArray, outArray)

  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface,   ONLY : Grid_getCellMetrics
  use Grid_data,        ONLY : gr_iMetricsGlb, gr_jMetricsGlb, gr_kMetricsGlb
  use gr_pfftInterface, ONLY : gr_pfftDcftForward, gr_pfftDcftInverse, gr_pfftTranspose, &
                               gr_pfftGetLocalLimitsAnytime, gr_pfftTriDiag, gr_pfftCyclicTriDiag
  use gr_pfftData,      ONLY : pfft_inLen, pfft_midLen, pfft_outLen, pfft_globalLen, pfft_transformType, &
                               pfft_work1, pfft_work2, pfft_procGrid, pfft_comm, pfft_myPE, pfft_scale,  &
                               pfft_trigIaxis, pfft_trigJaxis, pfft_trigKaxis

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Pfft.h"
#include "Flash_mpi.h"

  integer, intent(in)    :: iDirection, solveflag, inSize
  integer, dimension(MDIM), intent(IN) :: localSize, globalSize, transformType
  real, dimension(inSize), intent(IN)  :: inArray
  real, dimension(inSize), intent(OUT) :: outArray

  integer :: ierr
  integer :: J, M, N, size, nl
  integer, save :: ldw, liw, ilf, iuf
  integer, dimension(2,MDIM) :: pfftBlkLimits
  integer, dimension(:), allocatable, save :: iw
  real, save :: ch
  real, dimension(:), allocatable, save :: dw
  real, dimension(:), allocatable, save :: AM, BM, CM
  real, dimension(:), allocatable, save :: AN, BN, CN
  real, dimension(:,:), allocatable   :: temp2DArray, RHS
  real, dimension(:,:,:), allocatable :: temp3DArray
  logical, save :: firstCall = .true.
  logical, dimension(3), save :: init
  real :: mean, meanAux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                 !
!                                CONSTRUCT COEFFICIENT MATRICIES                                  !  
!                                                                                                 ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (firstCall) then


! Two Dimensional Solver
#if NDIM == 2

    ! --------------------------------------------------------------------------------------------------------!
    ! Initialize 2d block tridiagonal (pdc2d) ----------------------------------------------------------------!
    ! --------------------------------------------------------------------------------------------------------!

    ! dimensions
    N = globalSize(IAXIS) !NX-2  ! Total Number of Points in X 
    M = globalSize(JAXIS) !NY-2  ! Total Number of Points in Y 

    ! allocate coefficient arrays
    allocate(AN(N), BN(N), CN(N))
    allocate(AM(M), BM(M), CM(M))

    ! matrix coefficients for y-axis
    AM(1:M) = -gr_jMetricsGlb(LEFT_EDGE,1:M,1)
    BM(1:M) =  gr_jMetricsGlb(LEFT_EDGE,1:M,1) + gr_jMetricsGlb(RIGHT_EDGE,1:M,1)
    CM(1:M) = 1.0 / gr_jMetricsGlb(CENTER,1:M,1)
   
    ! apply boundary conditions y-axis
    BM(1) = gr_jMetricsGlb(RIGHT_EDGE,1,1)
    BM(M) = gr_jMetricsGlb(LEFT_EDGE, M,1)

    ! matrix coefficients for x-axis
    AN(1:N) = -gr_iMetricsGlb(LEFT_EDGE,1:N,1)
    BN(1:N) =  gr_iMetricsGlb(LEFT_EDGE,1:N,1) + gr_iMetricsGlb(RIGHT_EDGE,1:N,1)
    CN(1:N) = 1.0 / gr_iMetricsGlb(CENTER,1:N,1)

    ! apply boundary conditions x-axis
    BN(1) = gr_iMetricsGlb(RIGHT_EDGE,1,1)
    BN(N) = gr_iMetricsGlb(LEFT_EDGE, N,1)
 
    ! poisson coefficient
    ch = 0.0

    ! initialize solver
    init(1) = .false.
    init(2) = .true.
    init(3) = .true.

    ! initialize parameters for pdc2d
    size = pfft_procGrid(JAXIS)
    nl = 1 + max(int(log10(real(M))/log10(4.0)), 0)
    ldw = 6*nl*min((M + 2*size - 1)/size, M) + max(9*M, 11*N)
    liw = 6*M + (4**nl - 1)/3 + 2*nl + int(log10(real(2*size))/log10(4.0)) + 7
    allocate(dw(ldw), iw(liw))

    call pdc2dn(M, N, RHS, N, ilf, iuf, AM, BM, CM, AN, BN, CN, ch, dw, ldw, iw, liw, pfft_comm(JAXIS), init, ierr)

    ! prepare solver
    init(:) = .false.

    ! --------------------------------------------------------------------------------------------------------!
    ! Complete 2d block tridiagonal (pdc2d) ------------------------------------------------------------------!
    ! --------------------------------------------------------------------------------------------------------!


! Three Dimensional Solver 
#else

    ! --------------------------------------------------------------------------------------------------------!
    ! Initialize 3d block pentadiagonal (pdc3d) --------------------------------------------------------------!
    ! --------------------------------------------------------------------------------------------------------!
    

    ! ////// NOT YET IMPLEMENTED ///////// !


    ! --------------------------------------------------------------------------------------------------------!
    ! Complete 3d block pentadiagonal (pdc3d) ----------------------------------------------------------------!
    ! --------------------------------------------------------------------------------------------------------!

#endif

    ! Do not repeat initial setup
    firstCall = .false.

  endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                 !
!          PERFORM FORWARD DIRECTION TRANSPOSE AND MATRIX BASED LINEAR / PLANER SOLUTION          !  
!                                                                                                 ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (iDirection .eq. PFFT_FORWARD) then


! Two Dimensional Solver
#if NDIM == 2

    ! --------------------------------------------------------------------------------------------------------!
    ! Solve 2d block tridiagonal (pdc2d) ---------------------------------------------------------------------!
    ! --------------------------------------------------------------------------------------------------------!

    ! dimensions
    N = globalSize(IAXIS) !NX-2  ! Total Number of Points in X 
    M = globalSize(JAXIS) !NY-2  ! Total Number of Points in Y 

    ! lets work with the data as a 2d array {x,y} vice a 1d vector {x*y}
    allocate(RHS(pfft_inLen(IAXIS), pfft_inLen(JAXIS)))
    RHS = reshape(inArray, pfft_inLen(1:2))

    ! create solution array
    do J=1, pfft_inLen(JAXIS)
      RHS(:,J) = -(1.0/gr_iMetricsGlb(CENTER,1:N,1))*(1.0/gr_jMetricsGlb(CENTER,ilf+J-1,1))*RHS(:,J)
    end do
    
    ! solve the system
    call pdc2dn(M, N, RHS, N, ilf, iuf, AM, BM, CM, AN, BN, CN, ch, dw, ldw, iw, liw, pfft_comm(JAXIS), init, ierr)

    ! remove mean from zeroth wave component 
    meanAux = sum(RHS) / (M * N)
    call MPI_ALLreduce(meanAux, mean, 1, FLASH_REAL, MPI_SUM, pfft_comm(IAXIS), ierr)
    RHS(:,:) = RHS(:,:) - mean

    ! put the solution into the output array
    outArray(1:product(pfft_inLen)) = reshape(RHS, (/product(pfft_inLen)/))

    ! --------------------------------------------------------------------------------------------------------!
    ! Complete 2d block tridiagonal (pdc2d) ------------------------------------------------------------------!
    ! --------------------------------------------------------------------------------------------------------!


! Three Dimensional Solver
#else

    ! --------------------------------------------------------------------------------------------------------!
    ! Solve 3d block pentadiagonal (pdc3d) -------------------------------------------------------------------!
    ! --------------------------------------------------------------------------------------------------------!

    ! ////// NOT YET IMPLEMENTED ///////// !

    ! --------------------------------------------------------------------------------------------------------!
    ! Complete 2d block pentadiagonal (pdc3d)-----------------------------------------------------------------!
    ! --------------------------------------------------------------------------------------------------------!


#endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                 !
!                         PERFORM INVERSE DIRECTION TRANSPOSE                                     !  
!                                                                                                 ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  else


! Two Dimensional Solver
#if NDIM == 2

    ! --------------------------------------------------------------------------------------------------------!
    ! Solve 2d block tridiagonal (pdc2d) ---------------------------------------------------------------------!
    ! --------------------------------------------------------------------------------------------------------!

    ! ////// NOTHING TO IMPLEMENT FOR THIS METHOD ///////// !

    ! --------------------------------------------------------------------------------------------------------!
    ! Complete 2d block tridiagonal (pdc2d) ------------------------------------------------------------------!
    ! --------------------------------------------------------------------------------------------------------!

! Three Dimensional Solver
#else

    ! --------------------------------------------------------------------------------------------------------!
    ! Solve 3d block pentadiagonal (pdc3d) -------------------------------------------------------------------!
    ! --------------------------------------------------------------------------------------------------------!

    ! ////// NOT YET IMPLEMENTED ///////// !

    ! --------------------------------------------------------------------------------------------------------!
    ! Complete 2d block pentadiagonal (pdc3d)-----------------------------------------------------------------!
    ! --------------------------------------------------------------------------------------------------------!

#endif


  endif
  


  return

end subroutine gr_pfftPoissonDirect

