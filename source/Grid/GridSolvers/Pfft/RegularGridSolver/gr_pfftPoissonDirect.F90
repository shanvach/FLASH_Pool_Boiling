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
  use Grid_data,        ONLY : gr_iMetricsGlb, gr_jMetricsGlb, gr_kMetricsGlb
  use gr_pfftData,      ONLY : pfft_inLen, pfft_midLen, pfft_work1, pfft_work2, pfft_procGrid, pfft_trigIaxis, &
                               pfft_globalLen, pfft_comm, pfft_transformType, pfft_scale, pfft_myPE
!  use fish

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
  integer :: M, N, MP, NP
  real, dimension(:), allocatable, save :: AM, BM, CM, AN, BN, CN
!  TYPE (fishworkspace), save :: W
  real, dimension(:), allocatable :: RHS 
  real, dimension(:,:), allocatable :: temp2DArray
  logical, save :: firstCall = .true.
  real :: mean



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                 !
!                                CONSTRUCT COEFFICIENT MATRICIES                                  !  
!                                                                                                 ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (firstCall) then


! Two Dimensional Solver
#if NDIM == 2

    ! dimensions
    M = globalSize(IAXIS) !NX-2  ! Total Number of Points in X
    N = globalSize(JAXIS) !NY-2  ! Total Number of Points in Y

    ! allocate coefficient arrays
    allocate(AM(M), BM(M), CM(M), AN(N), BN(N), CN(N))

    ! matrix coefficients for blktri solve -- IAXIS
    AM(1:M) = gr_iMetricsGlb(CENTER,1:M,1) * gr_iMetricsGlb(LEFT_EDGE, 1:M,1)
    CM(1:M) = gr_iMetricsGlb(CENTER,1:M,1) * gr_iMetricsGlb(RIGHT_EDGE,1:M,1)
    MP = 0
    if (transformType(IAXIS) == PFFT_COS_CC) then
      if (pfft_myPE == 0) write(*,*) '2d pfft solver using IAXIS neumann coefficients'
      MP = 1
      AM(1) = 0.
      CM(M) = 0.
    endif
    BM = - AM - CM

    ! matrix coefficients for blktri solve -- JAXIS
    AN(1:N) = gr_jMetricsGlb(CENTER,1:N,1) * gr_jMetricsGlb(LEFT_EDGE, 1:N,1)
    CN(1:N) = gr_jMetricsGlb(CENTER,1:N,1) * gr_jMetricsGlb(RIGHT_EDGE,1:N,1)
    NP = 0
    if (transformType(JAXIS) == PFFT_COS_CC) then
      if (pfft_myPE == 0) write(*,*) '2d pfft solver using JAXIS neumann coefficients'
      NP = 1
      AN(1) = 0.
      CN(M) = 0.
    endif
    BN = - AN - CN

   ! allocate RHS array and blktri workspace
 
    

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

! Three Dimensional Solver
#else

    !! implement 2d transform + TDMA
    !! implement 1d transform + BlkTri

#endif


  endif
  


  return

end subroutine gr_pfftPoissonDirect

