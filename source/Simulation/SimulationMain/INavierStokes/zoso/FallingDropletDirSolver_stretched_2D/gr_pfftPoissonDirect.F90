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
  use Grid_data,        ONLY : gr_iMetricsGlb, gr_jMetricsGlb, gr_kMetricsGlb, gr_meshComm
  use gr_pfftInterface, ONLY : gr_pfftDcftForward, gr_pfftDcftInverse, gr_pfftTranspose, &
                               gr_pfftGetLocalLimitsAnytime, gr_pfftTriDiag, gr_pfftCyclicTriDiag
  use gr_pfftData,      ONLY : pfft_inLen, pfft_midLen, pfft_outLen, pfft_globalLen, pfft_transformType, &
                               pfft_work1, pfft_work2, pfft_procGrid, pfft_comm, pfft_myPE, pfft_scale,  &
                               pfft_trigIaxis, pfft_trigJaxis, pfft_trigKaxis, pfft_usableProc, pfft_me

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Pfft.h"
#include "Flash_mpi.h"

  integer, intent(in)    :: iDirection, solveflag, inSize
  integer, dimension(MDIM), intent(IN) :: localSize, globalSize, transformType
  real, dimension(inSize), intent(IN)  :: inArray
  real, dimension(inSize), intent(OUT) :: outArray

  integer :: ierr, error, errorAux, errorErr, errorMsg
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

  integer :: debugLen, i, k, gg, ii, jj, kk
  integer, dimension (:), allocatable  :: debugMsg, debugBuf
  real, dimension (:), allocatable  :: debugMsgF, debugBufF
  real, dimension(:,:), allocatable :: debugMsgArray, debugBufArray

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
 
    if (pfft_myPE == 0) write(*,*) '2d pfft solver -- 2d block tridiagonal (pdc2d)'

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
    select case (transformType(JAXIS))
    case (PFFT_COS_CC)
      if (pfft_myPE == 0) write(*,*) '2d pfft solver using JAXIS neumann coefficents'
      BM(1) = gr_jMetricsGlb(RIGHT_EDGE,1,1)
      BM(M) = gr_jMetricsGlb(LEFT_EDGE, M,1)
    case default
      call Driver_abortFlash("Unsupported y-direction boundary condition for RegularGridSolver!")
    end select

    ! matrix coefficients for x-axis
    AN(1:N) = -gr_iMetricsGlb(LEFT_EDGE,1:N,1)
    BN(1:N) =  gr_iMetricsGlb(LEFT_EDGE,1:N,1) + gr_iMetricsGlb(RIGHT_EDGE,1:N,1)
    CN(1:N) = 1.0 / gr_iMetricsGlb(CENTER,1:N,1)

    ! apply boundary conditions x-axis
    select case (transformType(IAXIS))
    case (PFFT_COS_CC)
      if (pfft_myPE == 0) write(*,*) '2d pfft solver using IAXIS neumann coefficents'
      BN(1) = gr_iMetricsGlb(RIGHT_EDGE,1,1)
      BN(N) = gr_iMetricsGlb(LEFT_EDGE, N,1)
    case default
      call Driver_abortFlash("Unsupported x-direction boundary condition for RegularGridSolver!")
    end select
 
    ! poisson coefficient
    ch = 0.0

    ! initialize solver
    init(1) = .false.
    init(2) = .true.
    init(3) = .true.

    ! initialize parameters for pdc2d
    call MPI_COMM_SIZE(pfft_comm(JAXIS), size, ierr)
    nl = 1 + max(int(log10(real(M))/log10(4.0)), 0)
    ldw = 6*nl*min((M + 2*size - 1)/size, M) + max(9*M, 11*N)
    liw = 6*M + (4**nl - 1)/3 + 2*nl + int(log10(real(2*size))/log10(4.0)) + 7
    allocate(dw(ldw), iw(liw))

    call pdc2dn(M, N, RHS, N, ilf, iuf, AM, BM, CM, AN, BN, CN, ch, dw, ldw, iw, liw, pfft_comm(JAXIS), init, ierr)

    ! identify if pdc2dn throws an error code
    errorAux = 0
    errorErr = 0
    if ( ierr .ne. 0 ) then 
      errorAux = 1
      errorErr = ierr
    endif
    !call MPI_ALLreduce(errorAux, error, 1, FLASH_REAL, MPI_SUM, pfft_comm(IAXIS), ierr)
    !call MPI_ALLreduce(errorErr, errorMsg, 1, FLASH_REAL, MPI_MAX, pfft_comm(IAXIS), ierr)
    if (error >= 1 .and. pfft_myPE == 0) then
      write(*,*) "Largest pdc2Dn generated error code among communicators is ", errorMsg
      call Driver_abortFlash("Error in pdc2Dn; internal solver error generated!")
    endif

    ! identify if pdc2dn strip matches pfft stripe in JAXIS 
    errorAux = 0
    errorErr = iuf - ilf + 1
    if ( pfft_inLen(JAXIS) .ne. errorErr ) then 
      errorAux = 1
    endif
    !call MPI_ALLreduce(errorAux, error, 1, FLASH_REAL, MPI_SUM, pfft_comm(IAXIS), ierr)
    !call MPI_ALLreduce(errorErr, errorMsg, 1, FLASH_REAL, MPI_MIN, pfft_comm(IAXIS), ierr)
    if (error >= 1 .and. pfft_myPE == 0) then
      write(*,*) "The pdc2Dn generated JAXIS stripe size does not match pfft; pdc2Dn min is ", errorMsg
      call Driver_abortFlash("Error in pdc2Dn; internal solver error generated!")
    endif

    ! Assure that other processes do not move forward
    call MPI_BARRIER(pfft_comm(IAXIS), ierr)

    ! prepare solver
    init(:) = .false.


    ! #$&*^%(^)*#_)*$)(&)($&)$_$&#)($&)$&
    ! DEBUG DEBUG DEBUG
    !if (pfft_myPE == 0) then
    ! do i=1, M
    !    print *, i, BM(i)
    !  end do
    !endif

    ! DEBUG DEBUG DEBUG
    ! #$&*^%(^)*#_)*$)(&)($&)$_$&#)($&)$&

    ! --------------------------------------------------------------------------------------------------------!
    ! Complete 2d block tridiagonal (pdc2d) ------------------------------------------------------------------!
    ! --------------------------------------------------------------------------------------------------------!


! Three Dimensional Solver 
#else

    ! --------------------------------------------------------------------------------------------------------!
    ! Initialize 3d block pentadiagonal (pdc3d) --------------------------------------------------------------!
    ! --------------------------------------------------------------------------------------------------------!
    
    if (pfft_myPE == 0) write(*,*) '3d pfft solver -- 3d block tridiagonal (pdc3d)'

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


    !######^#^#^#^#*#^*^#*^$*^&$(^&W$*(^&@$)(@&$)(@$&)($^&
    ! DEBUG  DEBUG  DEBUG

    debugLen = 1+8+3+4    
    call MPI_COMM_SIZE(pfft_comm(IAXIS), size, ierr)
    allocate(debugMsg(debugLen))
    if (pfft_myPE == 0) allocate(debugBuf(debugLen*size))    

    call gr_pfftGetLocalLimitsAnytime(IAXIS, IAXIS, IAXIS, pfft_inLen, PFFT_PCLDATA_REAL, pfftBlkLimits)
    call gr_pfftGetLocalLimitsAnytime(JAXIS, JAXIS, JAXIS, pfft_inLen, PFFT_PCLDATA_REAL, pfftBlkLimits)
    call gr_pfftGetLocalLimitsAnytime(KAXIS, KAXIS, KAXIS, pfft_inLen, PFFT_PCLDATA_REAL, pfftBlkLimits)
    
    call MPI_COMM_RANK(gr_meshComm, GG, ierr)
    call MPI_COMM_RANK(pfft_comm(IAXIS), II, ierr)
    call MPI_COMM_RANK(pfft_comm(JAXIS), JJ, ierr)
    call MPI_COMM_RANK(pfft_comm(KAXIS), KK, ierr)

!    write(*,*) GG
    i = 0
    if (pfft_usableProc) i = 1

    ! Gather debug info from processes
    debugMsg = (/ pfft_myPE,                                                &
                  GG, II, JJ, KK, i, pfft_me(1), pfft_me(2), pfft_me(3),        &
                  ilf, iuf, iuf-ilf+1,                                      &
                  pfftBlkLimits(LOW,JAXIS)+1,pfftBlkLimits(HIGH,JAXIS)+1,   &
                      pfftBlkLimits(HIGH,JAXIS)-pfftBlkLimits(LOW,JAXIS)+1, &
                      pfft_inLen(JAXIS)                                     &    
                /)
    call MPI_GATHER(debugMsg, debugLen, FLASH_INTEGER, debugBuf, debugLen, FLASH_INTEGER, 0, pfft_comm(IAXIS), ierr)
    
    !allocate(debugMsgArray(pfft_inLen(JAXIS),3)
    !if (pfft_myPE == 0) allocate(debugBufArray(pfft_inLen(JAXIS)*size,3))    
    !debugMsgArray(1:pfft_inLen(JAXIS),1) = AM(ilf:iuf)
    !debugMsgArray(1:pfft_inLen(JAXIS),2) = BM(ilf:iuf)
    !debugMsgArray(1:pfft_inLen(JAXIS),3) = CM(ilf:iuf)
    !call MPI_GATHER(debugMsgArray(:,1), pfft_inLen(JAXIS), FLASH_REAL, debugBufArray(:,1), pfft_inLen(JAXIS), FLASH_REAL, 0, pfft_comm(IAXIS), ierr)
    !call MPI_GATHER(debugMsgArray(:,2), pfft_inLen(JAXIS), FLASH_REAL, debugBufArray(:,2), pfft_inLen(JAXIS), FLASH_REAL, 0, pfft_comm(IAXIS), ierr)
    !call MPI_GATHER(debugMsgArray(:,3), pfft_inLen(JAXIS), FLASH_REAL, debugBufArray(:,3), pfft_inLen(JAXIS), FLASH_REAL, 0, pfft_comm(IAXIS), ierr)

    
    allocate(RHS(pfft_inLen(IAXIS), pfft_inLen(JAXIS)))
    RHS = reshape(inArray, pfft_inLen(1:2))
    allocate(debugMsgArray(pfft_inLen(JAXIS),pfft_inLen(IAXIS)))
    if (pfft_myPE == 0) allocate(debugBufArray(pfft_inLen(JAXIS)*size,pfft_inLen(IAXIS)))    
    do i=1, pfft_inLen(IAXIS)
      debugMsgArray(1:pfft_inLen(JAXIS),i) = RHS(i,1:pfft_inLen(JAXIS))           
      call MPI_GATHER(debugMsgArray(:,i), pfft_inLen(JAXIS), FLASH_REAL, debugBufArray(:,i), pfft_inLen(JAXIS), FLASH_REAL, 0, pfft_comm(IAXIS), ierr)
    end do
    deallocate(RHS)


    ! write debug info to stdout
    if ( pfft_myPE == 0 ) then
      print *, "N / M are", N, M 

      do i=0, size-1
        print '("0) Block ",I2," GRID ",I3,", IAXIS ",I3,", JAXIS ",I3,", KAXIS ",I3,", Use ",L," IAXIS ",I3,", JAXIS ",I3,", KAXIS ",I3)', &
          debugBuf(i*debugLen+1), debugBuf(i*debugLen+2:i*debugLen+9)
      end do
      print *, ""
      
      print *, "Block layout from pdc2d"
      do i=0, size-1    
        print '("1) Block ",I2," Range for JAXIS: ",I4," /",I4,"=",I4)', & 
          debugBuf(i*debugLen+1), debugBuf(i*debugLen+10:i*debugLen+12)
      end do
      print *, ""

      print *, "Block layout from pfft"
      do i=0, size-1    
        print '("2) Block ",I2," Range for JAXIS: ",I4," /",I4,"=",I4," /",I4)', &
          debugBuf(i*debugLen+1), debugBuf(i*debugLen+13:i*debugLen+16)
      end do
      print *, ""

      !ii = int(pfft_inLen(IAXIS) / 16)
      !jj = int(pfft_inLen(JAXIS) / 16)
      !do i=0, size-1
      !  print '("3) Block ",I4)', i
      !  do k=1, pfft_inLen(JAXIS), jj
      !     print '(*(f4.0,:,", "))', debugBufArray(i*pfft_inLen(JAXIS)+k,::ii)
      !  end do
      !end do 

    endif
    call MPI_BARRIER(pfft_comm(IAXIS), ierr)
    
    deallocate(debugMsg, debugMsgArray)
    if (pfft_myPE == 0) deallocate(debugBuf, debugBufArray)    
    ! DEBUG  DEBUG  DEBUG
    !######^#^#^#^#*#^*^#*^$*^&$(^&W$*(^&@$)(@&$)(@$&)($^&

    ! lets work with the data as a 2d array {x,y} vice a 1d vector {x*y}
    allocate(RHS(pfft_inLen(IAXIS), pfft_inLen(JAXIS)))
    RHS = reshape(inArray, pfft_inLen(1:2))
   
    ! create solution array
    do J=1, pfft_inLen(JAXIS) 
      RHS(:,J) = -(1.0/gr_iMetricsGlb(CENTER,1:N,1))*(1.0/gr_jMetricsGlb(CENTER,ilf+J-1,1))*RHS(:,J)
    end do
    !do J=pfftBlkLimits(LOW,JAXIS)+1, pfftBlkLimits(HIGH,JAXIS)+1 
    !  i = pfftBlkLimits(LOW,JAXIS)
    !  RHS(:,J-i) = -(1.0/gr_iMetricsGlb(CENTER,1:N,1))*(1.0/gr_jMetricsGlb(CENTER,J,1))*RHS(:,J-i)
    !end do
    
    ! solve the system
    call pdc2dn(M, N, RHS, N, ilf, iuf, AM, BM, CM, AN, BN, CN, ch, dw, ldw, iw, liw, pfft_comm(JAXIS), init, ierr)

    ! remove mean from zeroth wave component 
    meanAux = sum(RHS) / (M * N)
    call MPI_ALLreduce(meanAux, mean, 1, FLASH_REAL, MPI_SUM, pfft_comm(IAXIS), ierr)
    RHS(:,:) = RHS(:,:) - mean

    ! put the solution into the output array
    outArray(1:product(pfft_inLen)) = reshape(RHS, (/product(pfft_inLen)/))

    !######^#^#^#^#*#^*^#*^$*^&$(^&W$*(^&@$)(@&$)(@$&)($^&
    ! DEBUG  DEBUG  DEBUG
    
    debugLen = 1+2
    allocate(debugMsgF(debugLen))
    if (pfft_myPE == 0) allocate(debugBufF(debugLen*size))    

    debugMsgF = (/ float(pfft_myPE), minval(RHS), maxval(RHS) /)
    call MPI_GATHER(debugMsgF, debugLen, FLASH_REAL, debugBufF, debugLen, FLASH_REAL, 0, pfft_comm(IAXIS), ierr)

    ! write debug info to stdout
    if ( pfft_myPE == 0 ) then
      print '("MIN ",F6.2,", MAX ",F6.2)', minval(debugBufF(2::debugLen)), maxval(debugBufF(3::debugLen)) 

      do i=0, size-1
        print '("4) Block ",I3," MIN ",F6.2,", Max ",F6.2)', &
          int(debugBufF(i*debugLen+1)), debugBufF(i*debugLen+2:i*debugLen+3)
      end do
      print *, ""
    endif 

    deallocate(debugMsgF)
    if (pfft_myPE == 0) deallocate(debugBufF)    
    ! DEBUG  DEBUG  DEBUG
    !######^#^#^#^#*#^*^#*^$*^&$(^&W$*(^&@$)(@&$)(@$&)($^&

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

