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
  use gr_pfftData,      ONLY : pfft_inLen, pfft_midLen, pfft_outLen, pfft_globalLen, pfft_transformType, &
                               pfft_work1, pfft_work2, pfft_procGrid, pfft_comm, pfft_myPE, pfft_scale,  & 
                               pfft_trigIaxis, pfft_trigJaxis, pfft_trigKaxis,                           &
                               pfft_pclBaseDatType, pfft_t1Len, pfft_t2Len

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Pfft.h"
#include "Flash_mpi.h"

  integer, intent(in)    :: iDirection, solveflag, inSize
  integer, dimension(MDIM), intent(IN) :: localSize, globalSize, transformType
  real, dimension(inSize), intent(IN)  :: inArray
  real, dimension(inSize), intent(OUT) :: outArray

  integer :: numVec, ierr
  integer :: I, IL, J, JL, K, KL, ML, L, LL, M, N, NL, size
  integer, save :: ldw, liw, ilf, iuf
  integer, dimension(2,MDIM) :: pfftBlkLimits
  integer, dimension(:), allocatable, save :: iw
  real :: ch
  real, dimension(:), allocatable, save :: dw
  real, dimension(:), allocatable, save :: AM, BM, CM, CMM
  real, dimension(:), allocatable, save :: AN, BN, CN
  real, dimension(:), allocatable, save :: AK, AL, AJ
  real, dimension(:), allocatable :: BML, BMM, RHS, X 
  real, dimension(:,:), allocatable   :: temp2DArray
  real, dimension(:,:,:), allocatable :: temp3DArray
  logical, save :: firstCall = .true.
  logical, dimension(3), save :: init
  real :: mean, meanAux
  integer :: II, JJ, KK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                 !
!          CONSTRUCT COEFFICIENT MATRICIES AND DIAGONAL WAVE NUMBER COEFFICIENTS                  !  
!                                                                                                 ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (firstCall) then


! Two Dimensional Solver
#if NDIM == 2

    ! --------------------------------------------------------------------------------------------------------!
    ! Initialize 1d transform and 1d tridiagonal -------------------------------------------------------------!
    ! --------------------------------------------------------------------------------------------------------!
    
    ! dimensions
    L = globalSize(IAXIS) !NX-2  ! Total Number of Points in X
    M = globalSize(JAXIS) !NY-2  ! Total Number of Points in Y

    ! allocate coefficient arrays
    allocate(AK(L))
    allocate(AM(M), BM(M), CM(M))

    ! transform coefficients for x-axis
    select case (transformType(IAXIS)) 
    case (PFFT_REAL) 
      if (pfft_myPE == 0) write(*,*) '2d pfft solver using IAXIS periodic coefficents'
      do K=1, L/2
        AK(K) =  2.0 * PI * real(K-1)
      end do
      do K=L/2+1, L
        AK(K) = -2.0 * PI * real(L-K+1)
      end do
      AK(1:L) =  2.0 * ( 1.0 - cos(AK(1:L) / REAL(L)) ) * gr_iMetricsGlb(CENTER,1:L,1)**2
    case (PFFT_COS_CC)
      if (pfft_myPE == 0) write(*,*) '2d pfft solver using IAXIS neumann coefficents'
      AK(1) = 0.0
      do K=1, L-1
        AK(K+1) = 2.0 * ( 1.0 - cos(PI*real(K)/real(L)) ) * gr_iMetricsGlb(CENTER,K+1,1)**2
      end do
    case (PFFT_SIN_CC)
      if (pfft_myPE == 0) write(*,*) '2d pfft solver using IAXIS dirichlet coefficents'
      do K=1, L
        AK(K) = 2.0 * ( 1.0 - cos(PI*real(K)/real(L)) ) * gr_iMetricsGlb(CENTER,K,1)**2
      end do
    case default 
      call Driver_abortFlash("Unsupported x-direction transformation for RegularGridSolver!")
    end select

    ! matrix coefficients for tdma solve for y-axis
    AM(1:M) = gr_jMetricsGlb(CENTER,1:M,1) * gr_jMetricsGlb(LEFT_EDGE, 1:M,1)
    CM(1:M) = gr_jMetricsGlb(CENTER,1:M,1) * gr_jMetricsGlb(RIGHT_EDGE,1:M,1)
    select case (transformType(JAXIS))
    case (PFFT_REAL)
      if (pfft_myPE == 0) write(*,*) '2d pfft solver using JAXIS periodic coefficients'
    case (PFFT_COS_CC)
      if (pfft_myPE == 0) write(*,*) '2d pfft solver using JAXIS neumann coefficients'
      AM(1) = 0.0
      CM(M) = 0.0
    case (PFFT_SIN_CC)
      if (pfft_myPE == 0) write(*,*) '2d pfft solver using JAXIS dirichlet coefficients'
      AM(1) = 2.0 * AM(1)
      CM(M) = 2.0 * CM(M)
    case default
      call Driver_abortFlash("Unsupported y-direction transformation for RegularGridSolver!")
    end select
    BM = - AM - CM

    ! --------------------------------------------------------------------------------------------------------!
    ! Complete 1d transform and 1d tridiagonal----------------------------------------------------------------!
    ! --------------------------------------------------------------------------------------------------------!


! Three Dimensional Solver 
#else

    ! --------------------------------------------------------------------------------------------------------!
    ! initialize 2d transform and 1d tridiagonal !------------------------------------------------------------!
    if (solveflag .eq. TRIG_TRIG_DRCT) then      !------------------------------------------------------------!
      
      ! dimensions
      L = globalSize(IAXIS) !NX-1  ! Total Number of Points in X
      N = globalSize(JAXIS) !NY-1  ! Total Number of Points in Y
      M = globalSize(KAXIS) !NZ-1  ! Total Number of Points in Z

      ! allocate coefficient arrays
      allocate(AK(L))
      allocate(AL(N))
      allocate(AM(M), BM(M), CM(M))
    
      ! transformation Coefficients for the x-axis
      select case (transformType(IAXIS))
      case (PFFT_REAL)
        if (pfft_myPE == 0) write(*,*) '3d pfft solver using IAXIS periodic coefficents'
        do K=1, L/2
          AK(K) =  2.0 * PI * real(K-1)
        end do
        do K=L/2+1, L
          AK(K) = -2.0 * PI * real(L-K+1)
        end do
        AK(1:L) = 2.0 * ( 1.0 - cos(AK(1:L) / REAL(L)) ) * gr_iMetricsGlb(CENTER,1:L,1)**2
      case (PFFT_COS_CC)
        if (pfft_myPE == 0) write(*,*) '3d pfft solver using IAXIS neumann coefficents'
        AK(1) = 0.0
        do K=1, L-1
          AK(K+1) = 2.0 * ( 1.0 - cos(PI * real(K) / real(L)) ) * gr_iMetricsGlb(CENTER,K+1,1)**2
        end do
      case (PFFT_SIN_CC)
        if (pfft_myPE == 0) write(*,*) '2d pfft solver using IAXIS dirichlet coefficents'
        do K=1, L
          AK(K) = 2.0 * ( 1.0 - cos(PI * real(K) / real(L)) ) * gr_iMetricsGlb(CENTER,K,1)**2
        end do
      case default
        call Driver_abortFlash("Unsupported x-direction transformation for RegularGridSolver!")
      end select 

      ! transform coefficients for y-axis 
      select case (transformType(JAXIS))
      case (PFFT_REAL)
        if (pfft_myPE == 0) write(*,*) '3d pfft solver using JAXIS periodic coefficents'
        do K=1, N/2
          AL(K) =  2.0 * PI * real(K-1)
        end do
        do K=N/2+1, N
          AL(K) = -2.0 * PI * real(N-K+1)
        end do
        AL(1:N) = 2.0 * ( 1.0 - cos(AL(1:N) / REAL(N)) ) * gr_jMetricsGlb(CENTER,1:N,1)**2
      case (PFFT_COS_CC)
        if (pfft_myPE == 0) write(*,*) '3d pfft solver using JAXIS neumann coefficents'
        AL(1) = 0.0
        do K=1, N-1
          AL(K+1) = 2.0 * ( 1.0 - cos(PI * real(K) / real(N)) ) * gr_jMetricsGlb(CENTER,K+1,1)**2
        end do
      case (PFFT_SIN_CC)
        if (pfft_myPE == 0) write(*,*) '3d pfft solver using JAXIS dirichlet coefficents'
        do K=1, N
          AL(K) = 2.0 * ( 1.0 - cos(PI * real(K) / real(N)) ) * gr_jMetricsGlb(CENTER,K,1)**2
        end do
      case default 
        call Driver_abortFlash("Unsupported y-direction transformation for RegularGridSolver!")
      end select

      ! matrix coefficients for tdma solve for z-axis
      AM(1:M) = gr_kMetricsGlb(CENTER,1:M,1) * gr_kMetricsGlb(LEFT_EDGE, 1:M,1)
      CM(1:M) = gr_kMetricsGlb(CENTER,1:M,1) * gr_kMetricsGlb(RIGHT_EDGE,1:M,1)
      select case (transformType(KAXIS))
      case (PFFT_REAL)
        if (pfft_myPE == 0) write(*,*) '3d pfft solver using KAXIS periodic coefficents'
      case (PFFT_COS_CC)
        if (pfft_myPE == 0) write(*,*) '3d pfft solver using KAXIS neumann coefficents'
        AM(1) = 0.0
        CM(M) = 0.0
      case (PFFT_SIN_CC)
        if (pfft_myPE == 0) write(*,*) '3d pfft solver using KAXIS dirichlet coefficients'
        AM(1) = 2.0 * AM(1)
        CM(M) = 2.0 * CM(M)
      case default
        call Driver_abortFlash("Unsupported z-direction transformation for RegularGridSolver!")
      end select
      BM = - AM - CM
    
    ! --------------------------------------------------------------------------------------------------------!
    ! Complete 2d transform and 1d tridiagonal----------------------------------------------------------------!
    ! --------------------------------------------------------------------------------------------------------!


    ! --------------------------------------------------------------------------------------------------------!
    ! initialize 1d transform and 2d blktri (pdc2d)   !-------------------------------------------------------!
    else if (solveflag .eq. TRIG_DRCT_DRCT) then      !-------------------------------------------------------!

      ! dimensions
      L = globalSize(IAXIS) !NX-2  ! Total Number of Points in X
      N = globalSize(JAXIS) !NY-2  ! Total Number of Points in Y
      M = globalSize(KAXIS) !NZ-2  ! Total Number of Points in Z

      ! allocate coefficient arrays
      allocate(AK(L))
      allocate(AN(N), BN(N), CN(N))
      allocate(AM(M), BM(M), CM(M))

    ! transform coefficients for x-axis
    select case (transformType(IAXIS))
    case (PFFT_REAL)
      if (pfft_myPE == 0) write(*,*) '3d pfft solver using IAXIS periodic coefficents'
      do K=1, L/2
        AK(K) =  2.0 * PI * real(K-1)
      end do
      do K=L/2+1, L
        AK(K) = -2.0 * PI * real(L-K+1)
      end do
      AK(1:L) =  2.0 * ( 1.0 - cos(AK(1:L) / REAL(L)) ) * gr_iMetricsGlb(CENTER,1:L,1)**2
    case (PFFT_COS_CC)
      if (pfft_myPE == 0) write(*,*) '3d pfft solver using IAXIS neumann coefficents'
      AK(1) = 0.0
      do K=1, L-1
        AK(K+1) = 2.0 * ( 1.0 - cos(PI * real(K)/real(L)) ) * gr_iMetricsGlb(CENTER,K+1,1)**2
      end do
    case (PFFT_SIN_CC)
      if (pfft_myPE == 0) write(*,*) '3d pfft solver using IAXIS dirichlet coefficents'
      do K=1, L
        AK(K) = 2.0 * ( 1.0 - cos(PI*real(K)/real(L)) ) * gr_iMetricsGlb(CENTER,K,1)**2
      end do
    case default
      call Driver_abortFlash("Unsupported x-direction transformation for RegularGridSolver!")
    end select

    ! matrix coefficients for z-axis
    AM(1:M) = -gr_kMetricsGlb(LEFT_EDGE,1:M,1)
    BM(1:M) =  gr_kMetricsGlb(LEFT_EDGE,1:M,1) + gr_kMetricsGlb(RIGHT_EDGE,1:M,1)
    CM(1:M) = 1.0 / gr_kMetricsGlb(CENTER,1:M,1)

    ! apply boundary conditions z-axis
    BM(1) = gr_kMetricsGlb(RIGHT_EDGE,1,1)
    BM(M) = gr_kMetricsGlb(LEFT_EDGE, M,1)

    ! matrix coefficients for y-axis
    AN(1:N) = -gr_jMetricsGlb(LEFT_EDGE,1:N,1)
    BN(1:N) =  gr_jMetricsGlb(LEFT_EDGE,1:N,1) + gr_jMetricsGlb(RIGHT_EDGE,1:N,1)
    CN(1:N) = 1.0 / gr_jMetricsGlb(CENTER,1:N,1)

    ! apply boundary conditions y-axis
    BN(1) = gr_jMetricsGlb(RIGHT_EDGE,1,1)
    BN(N) = gr_jMetricsGlb(LEFT_EDGE, N,1)

    ! poisson coefficient
    ch = 0.0

    ! initialize solver
    init(1) = .false.
    init(2) = .true.
    init(3) = .true.

    ! initialize parameters for pdc2d
    call MPI_COMM_SIZE(pfft_comm(KAXIS), size, ierr) 
    nl = 1 + max(int(log10(real(M))/log10(4.0)), 0)
    ldw = 6*nl*min((M + 2*size - 1)/size, M) + max(9*M, 11*N)
    liw = 6*M + (4**nl - 1)/3 + 2*nl + int(log10(real(2*size))/log10(4.0)) + 7
    allocate(dw(ldw), iw(liw))

    call pdc2d(M, N, temp2DArray, N, ilf, iuf, AM, BM, CM, AN, BN, CN, ch, dw, ldw, iw, liw, pfft_comm(KAXIS), init, ierr)
    !write(*,*) "on process ", pfft_myPE, " ilf/iuf", ilf," ", iuf," span", pfft_midLen(JAXIS), "size", size    

    ! prepare solver
    init(:) = .false.

    ! --------------------------------------------------------------------------------------------------------!
    ! Complete 1d transform and 2d blktri  -------------------------------------------------------------------!
    ! --------------------------------------------------------------------------------------------------------!

    else
      call Driver_abortFlash("Unsupported Poisson Solver Encountered!")
    endif

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

    ! --------------------------------------------------------------------------------------------------------!
    ! Solve 1d transform and 1d tridiagonal ------------------------------------------------------------------!
    ! --------------------------------------------------------------------------------------------------------!

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
      if (transformType(JAXIS) == PFFT_COS_CC .or. transformType(JAXIS) == PFFT_SIN_CC) then
        call gr_pfftTriDiag(AM, BML, CM, RHS, X, M)     
      else 
        call gr_pfftCyclicTriDiag(AM, BML, CM, CM(M), AM(1), RHS, X, M)
      endif
     
      ! remove mean from zeroth wave component 
      if (J == 1) then
        mean = sum(X) / M      
        X(:) = X(:) - mean
      endif     
 
      ! store the solution
      temp2DArray(:,JL) = X(:)

    end do

    ! put the solution back into the pfft working array
    pfft_work2(1:product(pfft_midLen)) = reshape(temp2DArray, (/product(pfft_midLen)/))

    deallocate(temp2DArray)
    deallocate(BML, RHS, X)    

    ! --------------------------------------------------------------------------------------------------------!
    ! Complete 1d transform and 1d tridiagonal ---------------------------------------------------------------!
    ! --------------------------------------------------------------------------------------------------------!

! Three Dimensional Solver
#else

    ! --------------------------------------------------------------------------------------------------------!
    ! Solve 2d transform and 1d tridiagonal !-----------------------------------------------------------------!
    if (solveflag .eq. TRIG_TRIG_DRCT) then !-----------------------------------------------------------------!
      
      ! dimensions
      L = globalSize(IAXIS) !NX-2  ! Total Number of Points in X  
      N = globalSize(JAXIS) !NY-2  ! Total Number of Points in Y  
      M = globalSize(KAXIS) !NZ-2  ! Total Number of Points in Z  

      ! number of parallel transforms of x-direction 
      numVec = pfft_inLen(JAXIS) * pfft_inLen(KAXIS)

      ! forward transform the x-direction
      call gr_pfftDcftForward(inArray, pfft_work1, pfft_trigIaxis, pfft_globalLen(IAXIS), &
                              pfft_inLen(IAXIS), numVec, pfft_transformType(IAXIS), pfft_scale(IAXIS))

      ! transpose the pencil grid about y-axis :: {x,y,z} -> {y,z,x}
      call gr_pfftTranspose(iDirection, pfft_pclBaseDatType(IAXIS), pfft_work1, pfft_work2, &
                            pfft_t1Len, pfft_midLen, pfft_procGrid(JAXIS), pfft_comm(JAXIS))

      ! number of parallel transforms of y-direction
      numVec = pfft_midLen(JAXIS) * pfft_midLen(KAXIS)

      ! forward transform the y-direction
      call gr_pfftDcftForward(pfft_work2, pfft_work1, pfft_trigJaxis, pfft_globalLen(JAXIS), &
                              pfft_midLen(IAXIS), numVec, pfft_transformType(JAXIS), pfft_scale(JAXIS))

      ! transpose the pencil grid about z-axis :: {y,z,x} -> {z,x,y}
      call gr_pfftTranspose(iDirection, pfft_pclBaseDatType(JAXIS), pfft_work1, pfft_work2, &
                            pfft_t2Len, pfft_outLen, pfft_procGrid(KAXIS), pfft_comm(KAXIS))

      ! lets work with the data as a 3d array {z,x,y} vice a 1d vector {z*x*y}
      allocate(temp3DArray(pfft_outLen(IAXIS), pfft_outLen(JAXIS), pfft_outLen(KAXIS))) 
      temp3DArray = reshape(pfft_work2, pfft_outLen)
      
      allocate(BML(M), BMM(M), RHS(M), X(M))
      
      call gr_pfftGetLocalLimitsAnytime(IAXIS, KAXIS, 1, pfft_outLen, PFFT_PCLDATA_REAL, pfftBlkLimits)
      call gr_pfftGetLocalLimitsAnytime(JAXIS, IAXIS, 2, pfft_outLen, PFFT_PCLDATA_REAL, pfftBlkLimits)
      call gr_pfftGetLocalLimitsAnytime(KAXIS, JAXIS, 3, pfft_outLen, PFFT_PCLDATA_REAL, pfftBlkLimits)
      
      ! for each y-direction wave number 
      NL = pfft_outLen(KAXIS)
      do JL = 1, NL

        ! identify the wave number or y-direction index
        J = pfftBlkLimits(LOW,KAXIS) + JL
        if (J-1 > pfftBlkLimits(HIGH,KAXIS)) cycle

        ! center diagonal for a specific wave number in y-direction  
        select case (transformType(JAXIS))
        case (PFFT_COS_CC, PFFT_SIN_CC)
          BML(1:M) = BM(1:M) - AL(J)
        case (PFFT_REAL)
          BML(1:M) = BM(1:M) - AL(J/2+1)
        end select 

        ! for each x-direction wave number 
        LL = pfft_outLen(JAXIS)
        do IL = 1, LL

          ! identify the wave number or x-direction index
          I = pfftBlkLimits(LOW,JAXIS) + IL
          if (I-1 > pfftBlkLimits(HIGH,JAXIS)) cycle

          ! center diagonal for a specific wave number in x-direction  
          select case (transformType(IAXIS))
          case (PFFT_COS_CC, PFFT_SIN_CC)
            BMM(1:M) = BML(1:M) - AK(I)
          case (PFFT_REAL)
            BMM(1:M) = BML(1:M) - AK(I/2+1)
          end select

          ! right-hand-side and solution arrays
          RHS(:) = temp3DArray(:,IL,JL)
          X(:) = 0.

          ! solve the system
          if (transformType(KAXIS) == PFFT_COS_CC .or. transformType(KAXIS) == PFFT_SIN_CC) then
            call gr_pfftTriDiag(AM, BMM, CM, RHS, X, M)
          else
            call gr_pfftCyclicTriDiag(AM, BMM, CM, CM(M), AM(1), RHS, X, M)
          endif

          ! store the solution
          temp3DArray(:,IL,JL) = X(:)

        end do
      enddo

      ! put the solution back into the pfft working array
      pfft_work2(1:product(pfft_outLen)) = reshape(temp3DArray, (/product(pfft_outLen)/))
      deallocate(temp3DArray)
      deallocate(BMM, BML, RHS, X)

    ! --------------------------------------------------------------------------------------------------------!
    ! Complete 2d transform and 1d tridiagonal  --------------------------------------------------------------!
    ! --------------------------------------------------------------------------------------------------------!


    ! --------------------------------------------------------------------------------------------------------!
    ! Solve 1d transform and 2d blktri           !------------------------------------------------------------!
    else if (solveflag .eq. TRIG_DRCT_DRCT) then !------------------------------------------------------------!

      ! dimensions
      L = globalSize(IAXIS) !NX-2  ! Total Number of Points in X
      M = globalSize(KAXIS) !NZ-2  ! Total Number of Points in Z
      N = globalSize(JAXIS) !NY-2  ! Total Number of Points in Y

      ! number of parallel transforms of x-direction
      numVec = pfft_inLen(JAXIS) * pfft_inLen(KAXIS)

      ! forward transform the x-direction
      call gr_pfftDcftForward(inArray, pfft_work1, pfft_trigIaxis, pfft_globalLen(IAXIS), &
                              pfft_inLen(IAXIS), numVec, pfft_transformType(IAXIS), pfft_scale(IAXIS))

      ! transpose the pencil grid about y-axis :: {x,y,z} -> {y,z,x}
      call gr_pfftTranspose(iDirection, PFFT_PCLDATA_REAL, pfft_work1, pfft_work2, &
                            pfft_inLen, pfft_midLen, pfft_procGrid(JAXIS), pfft_comm(JAXIS))


      ! lets work with the data as a 3d array {y,z,x} vice a 1d vector {y*z*x}
      allocate(temp3DArray(pfft_midLen(IAXIS), pfft_midLen(JAXIS), pfft_midLen(KAXIS)))
      temp3DArray = reshape(pfft_work2, pfft_midLen(1:3))

      call gr_pfftGetLocalLimitsAnytime(IAXIS, JAXIS, 1, pfft_midLen, PFFT_PCLDATA_REAL, pfftBlkLimits)
      call gr_pfftGetLocalLimitsAnytime(JAXIS, KAXIS, 2, pfft_midLen, PFFT_PCLDATA_REAL, pfftBlkLimits)
      call gr_pfftGetLocalLimitsAnytime(KAXIS, IAXIS, 3, pfft_midLen, PFFT_PCLDATA_REAL, pfftBlkLimits)

      ! need to allocate a 2d right hand side for solver
      allocate(temp2DArray(pfft_midLen(IAXIS), pfft_midLen(JAXIS)))
      temp2DArray(:,:) = 0.0 

    !if (pfft_myPE == 0) write(*,*) "on process ", pfft_myPE, " in      ",pfft_inLen(IAXIS),  pfft_inLen(JAXIS),  pfft_inLen(KAXIS)
    !if (pfft_myPE == 0) write(*,*) "on process ", pfft_myPE, " mid     ",pfft_midLen(IAXIS), pfft_midLen(JAXIS), pfft_midLen(KAXIS)
    !if (pfft_myPE == 0) write(*,*) "on process ", pfft_myPE, " out     ",pfft_outLen(IAXIS), pfft_outLen(JAXIS), pfft_outLen(KAXIS)
    !call  MPI_COMM_RANK(pfft_comm(IAXIS), II, ierr)
    !call  MPI_COMM_RANK(pfft_comm(JAXIS), JJ, ierr)
    !call  MPI_COMM_RANK(pfft_comm(KAXIS), KK, ierr)
    !write(*,*) "on process ", pfft_myPE, " i/j/k   ", II, JJ, KK, "shape of 3dA", shape(temp3DArray), temp3DArray(1,1,1)


      ! for each x-direction wave number lets solve the block tridiagonal system in {y,z}
      LL = pfft_midLen(KAXIS)
      do JL = 1, LL

        ! identify the wave number or x-direction index
        J = pfftBlkLimits(LOW,JAXIS) + JL
        if (J-1 > pfftBlkLimits(HIGH,JAXIS)) cycle
        
        ! poisson coefficient
        select case (transformType(IAXIS))
        case (PFFT_COS_CC, PFFT_SIN_CC)
          ch = AK(J) 
        case (PFFT_REAL)
          ch = AK(J/2+1) 
        end select

        !if(pfft_myPE == 3) write(*,*) JL, J, ch
        !if(JL == 1) write(*,*) "on process ", pfft_myPE, ch
        
        ! create solution array
        do K=1, pfft_midLen(JAXIS)
          temp2DArray(:,K) = -(1.0/gr_jMetricsGlb(CENTER,1:N,1))*(1.0/gr_kMetricsGlb(CENTER,ilf+K-1,1))*temp3DArray(:,K,JL)
        end do

        ! solve the system
        call pdc2d(M, N, temp2DArray, N, ilf, iuf, AM, BM, CM, AN, BN, CN, ch, dw, ldw, iw, liw, pfft_comm(KAXIS), init, ierr)

        ! store the solution
        temp3DArray(:,:,JL) = temp2DArray(:,:)

      end do

      ! put the solution back into the pfft working array
      pfft_work2(1:product(pfft_midLen)) = reshape(temp3DArray, (/product(pfft_midLen)/))

      deallocate(temp3DArray)
      deallocate(temp2DArray)

    ! --------------------------------------------------------------------------------------------------------!
    ! Complete 1d transform and 2d blktri  -------------------------------------------------------------------!
    ! --------------------------------------------------------------------------------------------------------!


    else
      call Driver_abortFlash("Unsupported Poisson Solver Encountered!")
    endif

#endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                 !
!                         PERFORM INVERSE DIRECTION TRANSFORM                                     !  
!                                                                                                 ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  else


! Two Dimensional Solver
#if NDIM == 2

    ! --------------------------------------------------------------------------------------------------------!
    ! Solve 1d transform and 1d tridiagonal ------------------------------------------------------------------!
    ! --------------------------------------------------------------------------------------------------------!

    ! number of parallel transforms of x-direction
    numVec = pfft_inLen(JAXIS)

    ! transpose the pencil grid about y-axis :: {y,x} -> {x,y}
    call gr_pfftTranspose(iDirection, PFFT_PCLDATA_REAL, pfft_work2, pfft_work1, &
                          pfft_midLen, pfft_inLen, pfft_procGrid(JAXIS), pfft_comm(JAXIS))

    ! inverse transform the x-direction
    call gr_pfftDcftInverse(pfft_work1, outArray, pfft_trigIaxis, pfft_globalLen(IAXIS), &
                            pfft_inLen(IAXIS), numVec, pfft_transformType(IAXIS), 1.0)

    ! --------------------------------------------------------------------------------------------------------!
    ! Complete 1d transform and 1d tridiagonal ---------------------------------------------------------------!
    ! --------------------------------------------------------------------------------------------------------!

! Three Dimensional Solver
#else

    ! --------------------------------------------------------------------------------------------------------!
    ! Solve 2d transform and 1d tridiagonal  !----------------------------------------------------------------! 
    if (solveflag .eq. TRIG_TRIG_DRCT) then  !----------------------------------------------------------------!

      ! dimensions
      L = globalSize(IAXIS) !NX-2  ! Total Number of Points in X  
      N = globalSize(JAXIS) !NY-2  ! Total Number of Points in Y  
      M = globalSize(KAXIS) !NZ-2  ! Total Number of Points in Z  

      ! transpose the pencil grid about z-axis :: {z,x,y} -> {y,z,x}
      call gr_pfftTranspose(iDirection, pfft_pclBaseDatType(JAXIS), pfft_work2, pfft_work1, &
                            pfft_outLen, pfft_t2Len, pfft_procGrid(KAXIS), pfft_comm(KAXIS))
      
      ! number of parallel transforms of y-direction
      numVec = pfft_midLen(JAXIS) * pfft_midLen(KAXIS)

      ! inverse transform the y-direction
      call gr_pfftDcftInverse(pfft_work1, pfft_work2, pfft_trigJaxis, pfft_globalLen(JAXIS), &
                              pfft_midLen(IAXIS), numVec, pfft_transformType(JAXIS), 1.0) 

      ! transpose the pencil grid about y-axis :: {y,z,x} -> {x,y,z}
      call gr_pfftTranspose(iDirection, pfft_pclBaseDatType(IAXIS), pfft_work2, pfft_work1, &
                            pfft_midLen, pfft_t1Len, pfft_procGrid(JAXIS), pfft_comm(JAXIS))

      ! lets work with the data as a 3d array {x,y,z} vice a 1d vector {x*y*z}
      allocate(temp3DArray(pfft_inLen(IAXIS), pfft_inLen(JAXIS), pfft_inLen(KAXIS)))
      temp3DArray = reshape(pfft_work1, pfft_inLen)
      
      ! calculate mean from zeroth wave component 
      meanAux = sum(temp3DArray(1,:,:)) / (M * N)
      call MPI_ALLreduce(meanAux, mean, 1, FLASH_REAL, MPI_SUM, pfft_comm(IAXIS), ierr)
      
      ! remove mean from zeroth wave component 
      temp3DArray(1,:,:) = temp3DArray(1,:,:) - mean
      pfft_work1(1:product(pfft_inLen)) = reshape(temp3DArray, (/product(pfft_inLen)/))
      deallocate(temp3DArray)

      ! number of parallel transforms of x-direction 
      numVec = pfft_inLen(JAXIS) * pfft_inLen(KAXIS)

      ! forward transform the x-direction
      call gr_pfftDcftInverse(pfft_work1, outArray, pfft_trigIaxis, pfft_globalLen(IAXIS), &
                              pfft_inLen(IAXIS), numVec, pfft_transformType(IAXIS), 1.0)

    ! --------------------------------------------------------------------------------------------------------!
    ! Complete 2d transform and 1d tridiagonal  --------------------------------------------------------------!
    ! --------------------------------------------------------------------------------------------------------!

    ! --------------------------------------------------------------------------------------------------------!
    ! Solve 1d transform and 2d blktri           !------------------------------------------------------------!
    else if (solveflag .eq. TRIG_DRCT_DRCT) then !------------------------------------------------------------!

     ! dimensions
     L = globalSize(IAXIS) !NX-2  ! Total Number of Points in X
     M = globalSize(KAXIS) !NZ-2  ! Total Number of Points in Z
     N = globalSize(JAXIS) !NY-2  ! Total Number of Points in Y

    ! number of parallel transforms of x-direction
    numVec = pfft_inLen(JAXIS) * pfft_inLen(KAXIS)

    ! transpose the pencil grid about y-axis :: {y,z,x} -> {x,y,z}
    call gr_pfftTranspose(iDirection, PFFT_PCLDATA_REAL, pfft_work2, pfft_work1, &
                          pfft_midLen, pfft_inLen, pfft_procGrid(JAXIS), pfft_comm(JAXIS))

    ! lets work with the data as a 3d array {x,y,z} vice a 1d vector {x*y*z}
    allocate(temp3DArray(pfft_inLen(IAXIS), pfft_inLen(JAXIS), pfft_inLen(KAXIS)))
    temp3DArray = reshape(pfft_work1, pfft_inLen)

    

    ! calculate mean from zeroth wave component 
    meanAux = sum(temp3DArray(1,:,:)) / (M * N)
    call MPI_ALLreduce(meanAux, mean, 1, FLASH_REAL, MPI_SUM, pfft_comm(IAXIS), ierr)

    ! remove mean from zeroth wave component 
    temp3DArray(1,:,:) = temp3DArray(1,:,:) - mean
    pfft_work1(1:product(pfft_inLen)) = reshape(temp3DArray, (/product(pfft_inLen)/))
    deallocate(temp3DArray)
    
    !write(*,*) "meanAux on", pfft_myPE, " is", meanAux     
    !write(*,*) "mean on", pfft_myPE, " is", mean     

    ! inverse transform the x-direction
    call gr_pfftDcftInverse(pfft_work1, outArray, pfft_trigIaxis, pfft_globalLen(IAXIS), &
                            pfft_inLen(IAXIS), numVec, pfft_transformType(IAXIS), 1.0)

    ! --------------------------------------------------------------------------------------------------------!
    ! Complete 1d transform and 2d blktri  -------------------------------------------------------------------!
    ! --------------------------------------------------------------------------------------------------------!

    else
      call Driver_abortFlash("Unsupported Poisson Solver Encountered!")
    endif

#endif

  endif

  return

end subroutine gr_pfftPoissonTrigDirect

