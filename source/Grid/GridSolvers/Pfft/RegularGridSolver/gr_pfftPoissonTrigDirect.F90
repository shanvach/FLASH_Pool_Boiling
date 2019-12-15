!!****if* source/Grid/GridSolvers/Pfft/RegularGridSolver/gr_pfftPoissonTrigDirect
!!
!! NAME
!!
!!  gr_pfftPoissonTrigDirect
!!
!!
!! SYNOPSIS
!!
!!  call gr_pfftPoissonTrigDirect(integer(IN)                :: iDirection,
!!                                integer(IN)                :: solveflag,
!!                                integer(IN)                :: inSize,
!!                                integer(IN)                :: localSize(MDIM),
!!                                integer(IN)                :: globalSize(MDIM),
!!                                integer(IN)                :: transformType(MDIM),
!!                                real(IN) ,dimension(inSize):: inArray(:),
!!                                real(OUT),dimension(inSize):: outArray(:))
!!
!! DESCRIPTION
!!
!!   Poisson solver routine.  This module implements the mixed fft
!!   and direct based method for periodic, neumann, and dirichlet problems
!!   on a regular grid. Mixed boundaries for a given axis are not supported.
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

subroutine gr_pfftPoissonTrigDirect (iDirection, solveflag, inSize, localSize, globalSize, transformType, inArray, outArray)

  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface,   ONLY : Grid_getCellMetrics
  use Grid_data,        ONLY : gr_iMetricsGlb, gr_jMetricsGlb, gr_kMetricsGlb
  use gr_pfftInterface, ONLY : gr_pfftDcftForward, gr_pfftDcftInverse, gr_pfftTranspose, &
                               gr_pfftGetLocalLimitsAnytime, gr_pfftTriDiag, gr_pfftCyclicTriDiag
  use gr_pfftData,      ONLY : pfft_inLen, pfft_midLen, pfft_work1, pfft_work2, pfft_procGrid, pfft_trigIaxis, &
                               pfft_globalLen, pfft_comm, pfft_transformType, pfft_scale, pfft_myPE

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Pfft.h"
#include "Flash_mpi.h"

  integer, intent(in)    :: iDirection, solveflag, inSize
  integer, dimension(MDIM), intent(IN) :: localSize, globalSize, transformType
  real, dimension(inSize), intent(IN)  :: inArray
  real, dimension(inSize), intent(OUT) :: outArray

  integer :: numVec
  integer :: J, JL, K, L, LL, M
  integer, dimension(2,MDIM) :: pfftBlkLimits
  real, dimension(:), allocatable, save :: AM, BM, CM
  real, dimension(:), allocatable, save :: AK
  real, dimension(:), allocatable :: BML, RHS, X 
  real, dimension(:,:), allocatable :: temp2DArray
  logical, save :: firstCall = .true.
  real :: mean



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                 !
!          CONSTRUCT COEFFICIENT MATRICIES AND DIAGONAL WAVE NUMBER COEFFICIENTS                  !  
!                                                                                                 ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (firstCall) then


! Two Dimensional Solver
#if NDIM == 2

    ! dimensions
    L = globalSize(IAXIS) !NX-2  ! Total Number of Points in X
    M = globalSize(JAXIS) !NY-2  ! Total Number of Points in Y

    ! allocate coefficient arrays
    allocate(AK(L))
    allocate(AM(M), BM(M), CM(M))

    ! transform coefficients for x-axis (periodic transform -- fourier)
    if (transformType(IAXIS) == PFFT_REAL) then
      if (pfft_myPE == 0) write(*,*) '2d pfft solver using IAXIS periodic coefficents'
      do K=1, L/2
        AK(K) =  2. * PI * real(K-1)
      end do
      do K=L/2+1, L
        AK(K) = -2. * PI * real(L-K+1)
      end do
      AK(1:L) = 2. * ( 1. - cos(AK(1:L) / REAL(L)) ) * gr_iMetricsGlb(CENTER,1:L,1)**2

    ! transform coefficients for x-axis (neumann transform -- cos)
    else if (transformType(IAXIS) == PFFT_COS_CC) then
      if (pfft_myPE == 0) write(*,*) '2d pfft solver using IAXIS neumann coefficents'
      do K=1, L
        AK(K) =  1. * PI * real(K-1)
      end do
      AK(1:L) = 2. * ( 1. - cos(AK(1:L) / REAL(L)) ) * gr_iMetricsGlb(CENTER,1:L,1)**2

    ! unsupported transform type (e.g., SIN not currently supported)
    else 
      call Driver_abortFlash("Unsupported x-direction transformation for RegularGridSolver!")

    endif

    ! matrix coefficients for tdma solve
    AM(1:M) = gr_jMetricsGlb(CENTER,1:M,1) * gr_jMetricsGlb(LEFT_EDGE, 1:M,1)
    CM(1:M) = gr_jMetricsGlb(CENTER,1:M,1) * gr_jMetricsGlb(RIGHT_EDGE,1:M,1)
    if (transformType(JAXIS) == PFFT_COS_CC) then
      if (pfft_myPE == 0) write(*,*) '2d pfft solver using JAXIS neumann coefficients'
      AM(1) = 0.
      CM(M) = 0.
    endif
    BM = - AM - CM


! Three Dimensional Solver 
#else

    !! implement 2d transform + TDMA
    !! implement 1d transform + BlkTri

#endif

    ! Do not repeat initial setup
    firstCall = .false.

  endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                 !
!          PERFORM FORWARD DIRECTION TRANSFORM AND MATRIX BASED LINEAR / PLANER SOLUTION          !  
!                                                                                                 ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (iDirection .eq. PFFT_FORWARD) then


! Two Dimensional Solver
#if NDIM == 2

    ! dimensions
    L = globalSize(IAXIS) !NX-2  ! Total Number of Points in X
    M = globalSize(JAXIS) !NY-2  ! Total Number of Points in Y

    ! number of parallel transforms of x-direction
    numVec = pfft_inLen(JAXIS)

    ! forward transform the x-direction
    call gr_pfftDcftForward(inArray, pfft_work1, pfft_trigIaxis, pfft_globalLen(IAXIS), &
                            pfft_inLen(IAXIS), numVec, pfft_transformType(IAXIS), pfft_scale(IAXIS))

    ! transpose the pencil grid about y-axis :: {x,y} -> {y,x}
    call gr_pfftTranspose(iDirection, PFFT_PCLDATA_REAL, pfft_work1, pfft_work2, &
                          pfft_inLen, pfft_midLen, pfft_procGrid(JAXIS), pfft_comm(JAXIS))
 
    ! lets work with the data as a 2d array {y,x} vice a 1d vector {y*x}
    allocate(temp2DArray(pfft_midLen(IAXIS), pfft_midLen(JAXIS)))
    temp2DArray = reshape(pfft_work2, pfft_midLen(1:2))

    allocate(BML(M), RHS(M), X(M))

    call gr_pfftGetLocalLimitsAnytime(IAXIS, JAXIS, 1, pfft_midLen, PFFT_PCLDATA_REAL, pfftBlkLimits) 
    call gr_pfftGetLocalLimitsAnytime(JAXIS, IAXIS, 2, pfft_midLen, PFFT_PCLDATA_REAL, pfftBlkLimits)
    call gr_pfftGetLocalLimitsAnytime(KAXIS, KAXIS, 3, pfft_midLen, PFFT_PCLDATA_REAL, pfftBlkLimits)

    ! for each x-direction wave number lets solve the tridiagonal system in y
    LL = pfft_midLen(JAXIS)
    do JL = 1, LL

      ! identify the wave number or x-direction index
      J = pfftBlkLimits(LOW,JAXIS) + JL
      if (J-1 > pfftBlkLimits(HIGH,JAXIS)) cycle

      ! center diagonal and rhs for a specific wave number in x-direction
      !                               (+ tridiag modifies upper diagonal)
      BML(1:M) = BM(1:M) - AK(J/2+1)
      RHS(:) = temp2DArray(:,JL)
      X(:) = 0.
    
      ! solve the system
      if (transformType(JAXIS) == PFFT_COS_CC) then
        call gr_pfftTriDiag(AM, BML, CM, RHS, X, M)     
      else 
        call gr_pfftCyclicTriDiag(AM, BML, CM, CM(M), AM(1), RHS, X, M)
      endif
     
      ! remove mean from zeroth wave component 
      if (j == 1) then
        mean = sum(X) / M      
        X(:) = X(:) - mean
      endif     
 
      ! store the solution
      temp2DArray(:,JL) = X(:)

    end do

    ! put the solution back into the pfft working array
    pfft_work2(1:product(pfft_midLen)) = reshape(temp2DArray, (/product(pfft_midLen)/))

    deallocate(temp2DArray)
    deallocate(BML)    
    deallocate(RHS)    
    deallocate(X)    


! Three Dimensional Solver
#else

    !! implement 2d transform + TDMA
    !! implement 1d transform + BlkTri

#endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                 !
!                         PERFORM INVERSE DIRECTION TRANSFORM                                     !  
!                                                                                                 ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  else


! Two Dimensional Solver
#if NDIM == 2

    ! number of parallel transforms of x-direction
    numVec = pfft_inLen(JAXIS)

    ! transpose the pencil grid about y-axis :: {y,x} -> {x,y}
    call gr_pfftTranspose(iDirection, PFFT_PCLDATA_REAL, pfft_work2, pfft_work1, &
                          pfft_midLen, pfft_inLen, pfft_procGrid(JAXIS), pfft_comm(JAXIS))

    ! inverse transform the x-direction
    call gr_pfftDcftInverse(pfft_work1, outArray, pfft_trigIaxis, pfft_globalLen(IAXIS), &
                            pfft_inLen(IAXIS), numVec, pfft_transformType(IAXIS), 1.0)


! Three Dimensional Solver
#else

    !! implement 2d transform + TDMA
    !! implement 1d transform + BlkTri

#endif


  endif
  


  return

end subroutine gr_pfftPoissonTrigDirect

