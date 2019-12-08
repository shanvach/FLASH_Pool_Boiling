


subroutine gr_pfftPoissonDirect (iDirection, solveflag, inSize, localSize, globalSize, transformType, inArray, outArray)

  use Grid_interface,   ONLY : Grid_getCellMetrics
  use Grid_data,        ONLY : gr_iMetricsGlb, gr_jMetricsGlb, gr_kMetricsGlb
  use gr_pfftInterface, ONLY : gr_pfftDcftForward, gr_pfftDcftInverse, gr_pfftTranspose, &
                               gr_pfftGetLocalLimitsAnytime, gr_pfftTriDiag, gr_pfftCyclicTriDiag
  use gr_pfftData,      ONLY : pfft_inLen, pfft_midLen, pfft_work1, pfft_work2, pfft_procGrid, pfft_trigIaxis, &
                               pfft_globalLen, pfft_comm, pfft_transformType, pfft_scale, pfft_myPE
  use gr_interface,     ONLY : gr_findMean                             

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

  if (firstCall) then

     ! dimensions
     L = globalSize(IAXIS) !NX-2  ! Total Number of Points in X
     M = globalSize(JAXIS) !NY-2  ! Total Number of Points in Y

     ! allocate coefficient arrays
     allocate(AK(L))
     allocate(AM(M), BM(M), CM(M))

     if (transformType(IAXIS) == PFFT_REAL) then

       if (pfft_myPE == 0) write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!use periodic coefficents'
       ! transform coefficients for x-axis (periodic transform -- fourier)
       do K=1, L/2
          AK(K) =  2. * PI * real(K-1)
       end do
       do K=L/2+1, L
          AK(K) = -2. * PI * real(L-K+1)
       end do
       AK(1:L) = 2. * ( 1. - cos(AK(1:L) / REAL(L)) ) * gr_iMetricsGlb(CENTER,1:L,1)**2

     else

       if (pfft_myPE == 0) write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!use neumann coefficents'
       ! transform coefficients for x-axis (neumann transform -- cos)
       do K=1, L
          AK(K) =  1. * PI * real(K-1)
       end do
       AK(1:L) = 2. * ( 1. - cos(AK(1:L) / REAL(L)) ) * gr_iMetricsGlb(CENTER,1:L,1)**2

     endif

     ! matrix coefficients for tdma
     AM(1:M) = gr_jMetricsGlb(CENTER,1:M,1) * gr_jMetricsGlb(LEFT_EDGE, 1:M,1)
     CM(1:M) = gr_jMetricsGlb(CENTER,1:M,1) * gr_jMetricsGlb(RIGHT_EDGE,1:M,1)
     if (.true.) then
       if (pfft_myPE == 0) write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!zero out AM(1) and CM(M)'
       !AM(1) = 0.
       !CM(M) = 0.
     endif
     BM = - AM - CM

     firstCall = .false.

  endif


  if (iDirection .eq. PFFT_FORWARD) then

    !if (pfft_myPE == 0) write(*,*) 'begin forward solve'

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
        !if (pfft_myPE == 0) write(*,*) 'solve for ', JL
        call gr_pfftTriDiag(AM, BML, CM, RHS, X, M)     
      else 
        !if (pfft_myPE == 0) write(*,*) 'solve - cycle for ', JL
        call gr_pfftCyclicTriDiag(AM, BML, CM, CM(M), AM(1), RHS, X, M)
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

  else


    !if (pfft_myPE == 0) write(*,*) 'begin inverse solve'

    ! number of parallel transforms of x-direction
    numVec = pfft_inLen(JAXIS)

    ! transpose the pencil grid about y-axis :: {y,x} -> {x,y}
    call gr_pfftTranspose(iDirection, PFFT_PCLDATA_REAL, pfft_work2, pfft_work1, &
                          pfft_midLen, pfft_inLen, pfft_procGrid(JAXIS), pfft_comm(JAXIS))

    ! inverse transform the x-direction
    call gr_pfftDcftInverse(pfft_work1, outArray, pfft_trigIaxis, pfft_globalLen(IAXIS), &
                            pfft_inLen(IAXIS), numVec, pfft_transformType(IAXIS), 1.0)

  endif
  
  return

end subroutine gr_pfftPoissonDirect

