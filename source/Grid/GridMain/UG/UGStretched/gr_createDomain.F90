!!****if* source/Grid/GridMain/UG/gr_createDomain
!!
!! NAME
!!
!!  gr_createDomain
!!
!!
!! SYNOPSIS
!!
!!  gr_createDomain()
!!                  
!!
!!
!! DESCRIPTION
!!
!! Creates the Uniform Grid domain. The Uniform Grid expects to be
!! given the processor grid at runtime in the form of runtime parameters
!! iprocs,jprocs and  kprocs, where iprocs*jprocs*kprocs is the
!! number of processor in the run. If UG is running in fixed
!! blocksize mode, then the blocksizes NXB, NYB and NZB are specified
!! at compile time. One block is placed on each processor, so that
!! the global domain size is <NXB*iprocs, NYB*jprocs, NZB*kprocs>.
!! However, if UG is running in non-fixed blocksize mode, the
!! blocksize is determined at runtime. In this mode, UG expects to be
!! given the global domain size in form of runtime parameters
!! iGridSize,jGridSize and kGridSize, and the blocksize is <iGridSize
!!/iprocs,jGridSize/jprocs,kGridSize/kprocs>. As in fixedblocksize mode,
!! only one blocks is placed on each processor
!!  
!! This routine also creates directional communicators used in
!! guardcell exchanges, allocates space for storing physical
!! coordinates of grid points, and computes them.
!!
!!  
!! ARGUMENTS
!!
!!
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine gr_createDomain()

  use physicalData, ONLY : unk, facevarx,facevary,facevarz

  use Grid_data, ONLY : scratch,scratch_ctr,&
       &scratch_facevarx,scratch_facevary,scratch_facevarz
  
  use Grid_data, ONLY : gr_axisMe, gr_axisNumProcs,&
       gr_guard, &
       gr_gIndexSize, gr_blkCornerID, &
       gr_lIndexSize, gr_kmax, gr_kmin, &
       gr_jmax, gr_jmin, gr_imax, gr_imin, gr_delta, gr_blkBC, &
       gr_iCoords, gr_jCoords, gr_kCoords, gr_domainBC, gr_iLoGC, &
       gr_jLoGC, gr_kLoGC, gr_iHiGC, gr_jHiGC, gr_kHiGC,&
       gr_ilo,gr_ihi,gr_jlo,gr_jhi,gr_klo,gr_khi,&
       gr_iguard,gr_jguard,gr_kguard,&
       gr_iStr,gr_jStr,gr_kStr,&
       gr_iStrType,gr_jStrType,gr_kStrType,&
       gr_iStrPar,gr_jStrPar,gr_kStrPar,&
       gr_iMetrics,gr_jMetrics,gr_kMetrics,&
       gr_iMetricsGlb,gr_jMetricsGlb,gr_kMetricsGlb,&
       gr_iCoordsGlb, gr_jCoordsGlb, gr_kCoordsGlb
  
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer :: ierr
  integer :: color, key, range_b, range_e
  real :: halfDelta
  integer :: i,j
  integer :: gI, gJ, gK

  !store local index size for each block
  gr_lIndexSize = gr_gIndexSize/gr_axisNumProcs

  !store lower left global index for each dim
  gr_blkCornerID = gr_axisMe*gr_lIndexSize+1  

  gr_iloGc = 1
  gr_ilo = gr_iloGc+gr_guard(IAXIS)
  gr_ihi = gr_ilo+gr_lIndexSize(IAXIS) -1
  gr_ihiGc = gr_ihi + gr_guard(IAXIS)
  gr_iguard = gr_guard(IAXIS)
  gI = gr_gIndexSize(IAXIS)

  if(NDIM>1)then

     gr_jloGc = 1
     gr_jlo = gr_jloGc+gr_guard(JAXIS)
     gr_jhi = gr_jlo+gr_lIndexSize(JAXIS) -1
     gr_jhiGc = gr_jhi + gr_guard(JAXIS)
     gr_jguard = gr_guard(JAXIS)
     gJ = gr_gIndexSize(JAXIS)

  else
     gr_jloGc = 1
     gr_jlo = 1
     gr_jhi = 1
     gr_jhiGc = 1
     gr_jguard = 0
  end if


  if(NDIM>2)then
     gr_kloGc = 1
     gr_klo = gr_kloGc+gr_guard(KAXIS)
     gr_khi = gr_klo+gr_lIndexSize(KAXIS) -1
     gr_khiGc = gr_khi + gr_guard(KAXIS)
     gr_kguard = gr_guard(KAXIS)
     gK = gr_gIndexSize(KAXIS)

  else
     gr_klo=1
     gr_kloGc=1
     gr_khi = 1
     gr_khiGc = 1
     gr_kguard = 0
  endif
  
  !! Now create the grid and coordinates etc
  allocate(gr_iCoords(3,gr_ihiGc-gr_iloGc+1,1))
  allocate(gr_jCoords(3,gr_jhiGc-gr_jloGc+1,1))
  allocate(gr_kCoords(3,gr_khiGc-gr_kloGc+1,1))

  allocate(gr_iMetrics(3,gr_ihiGc-gr_iloGc+1,1))
  allocate(gr_jMetrics(3,gr_jhiGc-gr_jloGc+1,1))
  allocate(gr_kMetrics(3,gr_khiGc-gr_kloGc+1,1))

  allocate(gr_iMetricsGlb(3,gI,1))
  allocate(gr_jMetricsGlb(3,gJ,1))
  allocate(gr_kMetricsGlb(3,gK,1))

  allocate(gr_iCoordsGlb(3,gI,1))
  allocate(gr_jCoordsGlb(3,gJ,1))
  allocate(gr_kCoordsGlb(3,gK,1))

#ifndef FIXEDBLOCKSIZE
  
  allocate(unk(UNK_VARS_BEGIN:UNK_VARS_END,&
       gr_iLoGc:gr_iHiGc, gr_jLoGc:gr_jHiGc,&
       gr_kLoGc:gr_kHiGc,1))
  
  
#if(NFACE_VARS>0)
  
  allocate(facevarx( NFACE_VARS,&
       gr_iLoGc:gr_iHiGc+1, gr_jLoGc:gr_jHiGc,&
       gr_kLoGc:gr_kHiGc,1))
  
  allocate(facevary( NFACE_VARS,&
       gr_iLoGc:gr_iHiGc, gr_jLoGc:gr_jHiGc+K2D,&
       gr_kLoGc:gr_kHiGc,1))
  
  allocate(facevarz( NFACE_VARS,&
       gr_iLoGc:gr_iHiGc, gr_jLoGc:gr_jHiGc,&
       gr_kLoGc:gr_kHiGc+K3D,1) )
  
#else
  allocate(facevarx(1,1,1,1,1))
  allocate(facevarz(1,1,1,1,1))
  allocate(facevary(1,1,1,1,1))
#endif
#if NSCRATCH_GRID_VARS > 0
  allocate(scratch(SCRATCH_GRID_VARS_BEGIN:SCRATCH_GRID_VARS_END,&
       gr_iLoGc:gr_iHiGc+1, gr_jLoGc:gr_jHiGc+1,&
       gr_kLoGc:gr_kHiGc+1,1))
#else
  allocate(scratch(1,1,1,1,1))
#endif

#if NSCRATCH_CENTER_VARS > 0
  allocate(scratch_ctr(SCRATCH_CENTER_VARS_BEGIN:SCRATCH_CENTER_VARS_END,&
       gr_iLoGc:gr_iHiGc, gr_jLoGc:gr_jHiGc,&
       gr_kLoGc:gr_kHiGc,1))
#else
  allocate(scratch_ctr(1,1,1,1,1))
#endif

#if(NSCRATCH_FACEX_VARS>0)  
  allocate(scratch_facevarx( SCRATCH_FACEX_VARS_BEGIN:SCRATCH_FACEX_VARS_END,&
       gr_iLoGc:gr_iHiGc+1, gr_jLoGc:gr_jHiGc,&
       gr_kLoGc:gr_kHiGc,1))
#else
  allocate(scratch_facevarx(1,1,1,1,1))
#endif

#if(NSCRATCH_FACEY_VARS>0)  
  allocate(scratch_facevary( SCRATCH_FACEY_VARS_BEGIN:SCRATCH_FACEY_VARS_END,&
       gr_iLoGc:gr_iHiGc, gr_jLoGc:gr_jHiGc+K2D,&
       gr_kLoGc:gr_kHiGc,1))
#else
  allocate(scratch_facevary(1,1,1,1,1))
#endif  

#if(NSCRATCH_FACEZ_VARS>0)
  allocate(scratch_facevarz( SCRATCH_FACEZ_VARS_BEGIN:SCRATCH_FACEZ_VARS_END,&
       gr_iLoGc:gr_iHiGc, gr_jLoGc:gr_jHiGc,&
       gr_kLoGc:gr_kHiGc+K3D,1) )
#else
  allocate(scratch_facevarz(1,1,1,1,1))
#endif
  
#endif

#ifdef DEBUG_GRID
  write(6,*)'gr_gIndexSize', gr_gIndexSize
  write(6,*)'gr_lIndexSize', gr_lIndexSize
  write(6,*)'gr_blkCornerID',gr_blkCornerID
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Create x grid coordinates and metrics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  if(gr_iStr) then !stretched grid in x
    
     ! grid stretching in x axis not currently supported
     print*,"Create Domain : invalid stretching axis "
     call Driver_abortFlash("Create Domain : invalid stretching axis ")

  else ! uniform grid

    ! calculate block cell IAXIS coordinate
    gr_delta(IAXIS) = (gr_imax-gr_imin)/gr_gIndexSize(IAXIS)
    halfDelta = gr_delta(IAXIS)/2.0
    j = gr_blkCornerID(IAXIS)-gr_guard(IAXIS)-1
    
    do i = gr_iloGc,gr_ihiGc
       gr_iCoords(LEFT_EDGE,i,1) = gr_imin+j*gr_delta(IAXIS)
       gr_iCoords(CENTER,i,1) = gr_imin+j*gr_delta(IAXIS)+halfDelta
       j = j+1
       gr_iCoords(RIGHT_EDGE,i,1) = gr_imin+j*gr_delta(IAXIS)
    end do

    ! calculate global cell IAXIS coordinates
    j = 0
    do i = 1,gI
       gr_iCoordsGlb(LEFT_EDGE,i,1) = gr_imin+j*gr_delta(IAXIS)
       gr_iCoordsGlb(CENTER,i,1) = gr_imin+j*gr_delta(IAXIS)+halfDelta
       j = j+1
       gr_iCoordsGlb(RIGHT_EDGE,i,1) = gr_imin+j*gr_delta(IAXIS)
    end do

  end if
 
  gr_blkBC = gr_domainBC
  if(gr_axisMe(IAXIS)/=0)gr_blkBC(LOW,IAXIS)=NOT_BOUNDARY
  if(gr_axisMe(IAXIS)/=(gr_axisNumProcs(IAXIS)-1))&
     gr_blkBC(HIGH,IAXIS)=NOT_BOUNDARY
  
  ! calculate block IAXIS direction metric coefficients
  gr_iMetrics(CENTER,:,1)                       = 1.0 / (gr_iCoords(RIGHT_EDGE,gr_iloGC:gr_ihiGC,1) - gr_iCoords(LEFT_EDGE,gr_iloGC:gr_ihiGC,1))
  gr_iMetrics(RIGHT_EDGE,gr_iloGC:gr_ihiGC-1,1) = 1.0 / (gr_iCoords(CENTER,gr_iloGC+1:gr_ihiGC,1)   - gr_iCoords(CENTER,gr_iloGC:gr_ihiGC-1,1))
  gr_iMetrics(RIGHT_EDGE,gr_ihiGC,1)            = 2.0 / (-3.0*gr_iCoords(RIGHT_EDGE,gr_ihiGC,1) + 4.0*gr_iCoords(RIGHT_EDGE,gr_ihiGC-1,1) - gr_iCoords(RIGHT_EDGE,gr_ihiGC-2,1))
  gr_iMetrics(LEFT_EDGE, gr_iloGC+1:gr_ihiGC,1) = 1.0 / (gr_iCoords(CENTER,gr_iloGC+1:gr_ihiGC,1)   - gr_iCoords(CENTER,gr_iloGC:gr_ihiGC-1,1))
  gr_iMetrics(LEFT_EDGE, gr_iloGC,1)            = 2.0 / (-3.0*gr_iCoords(LEFT_EDGE,gr_iloGC,1)  + 4.0*gr_iCoords(LEFT_EDGE,gr_iloGC+1,1)  - gr_iCoords(LEFT_EDGE,gr_iloGC+2,1))

  ! calculate global IAXIS direction metric coefficients
  gr_iMetricsGlb(CENTER,1:gI,1)       = 1.0 / (gr_iCoordsGlb(RIGHT_EDGE,1:gI,1) - gr_iCoordsGlb(LEFT_EDGE,1:gI,1))
  gr_iMetricsGlb(RIGHT_EDGE,1:gI-1,1) = 1.0 / (gr_iCoordsGlb(CENTER,2:gI,1)     - gr_iCoordsGlb(CENTER,1:gI-1,1))
  gr_iMetricsGlb(RIGHT_EDGE,gI,1)     = 2.0 / (-3.0*gr_iCoordsGlb(RIGHT_EDGE,gI,1) + 4.0*gr_iCoordsGlb(RIGHT_EDGE,gI-1,1) - gr_iCoordsGlb(RIGHT_EDGE,gI-2,1))
  gr_iMetricsGlb(LEFT_EDGE, 2:gI,1)   = 1.0 / (gr_iCoordsGlb(CENTER,2:gI,1)     - gr_iCoordsGlb(CENTER,1:gI-1,1))
  gr_iMetricsGlb(LEFT_EDGE, 1,1)      = 2.0 / (-3.0*gr_iCoordsGlb(LEFT_EDGE,1,1)   + 4.0*gr_iCoordsGlb(LEFT_EDGE,2,1)     - gr_iCoordsGlb(LEFT_EDGE,3,1))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Create y grid coordinates and metrics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(NDIM > 1) then ! 2d or 3d simulation

    if(gr_jStr) then !stretched grid in y
      select case (gr_jStrType)
      
      case(SG_USER)

        ! stretching method SG_USER not currently supported
        print*,"Create Domain : invalid stretching method "
        call Driver_abortFlash("Create Domain : invalid stretching method ")

      case(SG_TANH)

	! calculate block cell JAXIS coordinates
        ! start with a uniform distribution of points [0,1]
        gr_delta(JAXIS) = 1.0/gr_gIndexSize(JAXIS)
        halfDelta = gr_delta(JAXIS)/2.0
        j = gr_blkCornerID(JAXIS)-gr_guard(JAXIS)-1
        
        ! transform uniform to stretched using y(s) = Tanh((-1+2s) atan(a)) + 1 / 2a
        do i = gr_jloGc,gr_jhiGc
          gr_jCoords(LEFT_EDGE,i,1)  = (gr_jmax-gr_jmin)*(tanh((-1.0+2.0*j*gr_delta(JAXIS)          )*atanh(gr_jStrPar))+1.0)/(2.0*gr_jStrPar)+gr_jmin
          gr_jCoords(CENTER,i,1)     = (gr_jmax-gr_jmin)*(tanh((-1.0+2.0*j*gr_delta(JAXIS)+halfDelta)*atanh(gr_jStrPar))+1.0)/(2.0*gr_jStrPar)+gr_jmin
          j = j+1
          gr_jCoords(RIGHT_EDGE,i,1) = (gr_jmax-gr_jmin)*(tanh((-1.0+2.0*j*gr_delta(JAXIS)          )*atanh(gr_jStrPar))+1.0)/(2.0*gr_jStrPar)+gr_jmin
        end do

        ! calculate global cell JAXIS coordinates
        j = 0
        do i = 1,gJ
          gr_jCoordsGlb(LEFT_EDGE,i,1)  = (gr_jmax-gr_jmin)*(tanh((-1.0+2.0*j*gr_delta(JAXIS)          )*atanh(gr_jStrPar))+1.0)/(2.0*gr_jStrPar)+gr_jmin
          gr_jCoordsGlb(CENTER,i,1)     = (gr_jmax-gr_jmin)*(tanh((-1.0+2.0*j*gr_delta(JAXIS)+halfDelta)*atanh(gr_jStrPar))+1.0)/(2.0*gr_jStrPar)+gr_jmin
          j = j+1
          gr_jCoordsGlb(RIGHT_EDGE,i,1) = (gr_jmax-gr_jmin)*(tanh((-1.0+2.0*j*gr_delta(JAXIS)          )*atanh(gr_jStrPar))+1.0)/(2.0*gr_jStrPar)+gr_jmin
        end do
      
      case default

        ! stretching method unknown
        print*,"Create Domain : invalid stretching method "
        call Driver_abortFlash("Create Domain : invalid stretching method ")

      end select

    else ! uniform grid

      ! calculate block cell JAXIS coordinates
      gr_delta(JAXIS) = (gr_jmax-gr_jmin)/gr_gIndexSize(JAXIS)
      halfDelta = gr_delta(JAXIS)/2.0
      j = gr_blkCornerID(JAXIS)-gr_guard(JAXIS)-1
  
      do i = gr_jloGc,gr_jhiGc
         gr_jCoords(LEFT_EDGE,i,1) = gr_jmin+j*gr_delta(JAXIS)
         gr_jCoords(CENTER,i,1) = gr_jmin+j*gr_delta(JAXIS)+halfDelta
         j = j+1
         gr_jCoords(RIGHT_EDGE,i,1) = gr_jmin+j*gr_delta(JAXIS)
      end do

      ! calculate global cell JAXIS coordinates
      j = 0
      do i = 1,gJ
         gr_jCoordsGlb(LEFT_EDGE,i,1) = gr_jmin+j*gr_delta(JAXIS)
         gr_jCoordsGlb(CENTER,i,1) = gr_jmin+j*gr_delta(JAXIS)+halfDelta
         j = j+1
         gr_jCoordsGlb(RIGHT_EDGE,i,1) = gr_jmin+j*gr_delta(JAXIS)
      end do

    end if

    if(gr_axisMe(JAXIS)/=0)gr_blkBC(LOW,JAXIS)=NOT_BOUNDARY
    if(gr_axisMe(JAXIS)/=(gr_axisNumProcs(JAXIS)-1))&
       gr_blkBC(HIGH,JAXIS)=NOT_BOUNDARY

    ! calculate block JAXIS direction metric coefficients
    gr_jMetrics(CENTER,:,1)                       = 1.0 / (gr_jCoords(RIGHT_EDGE,gr_jloGC:gr_jhiGC,1) - gr_jCoords(LEFT_EDGE,gr_jloGC:gr_jhiGC,1))
    gr_jMetrics(RIGHT_EDGE,gr_jloGC:gr_jhiGC-1,1) = 1.0 / (gr_jCoords(CENTER,gr_jloGC+1:gr_jhiGC,1)   - gr_jCoords(CENTER,gr_jloGC:gr_jhiGC-1,1))
    gr_jMetrics(RIGHT_EDGE,gr_jhiGC,1)            = 2.0 / (-3.0*gr_jCoords(RIGHT_EDGE,gr_jhiGC,1) + 4.0*gr_jCoords(RIGHT_EDGE,gr_jhiGC-1,1) - gr_jCoords(RIGHT_EDGE,gr_jhiGC-2,1))
    gr_jMetrics(LEFT_EDGE, gr_jloGC+1:gr_jhiGC,1) = 1.0 / (gr_jCoords(CENTER,gr_jloGC+1:gr_jhiGC,1)   - gr_jCoords(CENTER,gr_jloGC:gr_jhiGC-1,1))
    gr_jMetrics(LEFT_EDGE, gr_jloGC,1)            = 2.0 / (-3.0*gr_jCoords(LEFT_EDGE,gr_jloGC,1)  + 4.0*gr_jCoords(LEFT_EDGE,gr_jloGC+1,1)  - gr_jCoords(LEFT_EDGE,gr_jloGC+2,1))

    ! calculate global JAXIS direction metric coefficients
    gr_jMetricsGlb(CENTER,1:gJ,1)       = 1.0 / (gr_jCoordsGlb(RIGHT_EDGE,1:gJ,1) - gr_jCoordsGlb(LEFT_EDGE,1:gJ,1))
    gr_jMetricsGlb(RIGHT_EDGE,1:gJ-1,1) = 1.0 / (gr_jCoordsGlb(CENTER,2:gJ,1)     - gr_jCoordsGlb(CENTER,1:gJ-1,1))
    gr_jMetricsGlb(RIGHT_EDGE,gJ,1)     = 2.0 / (-3.0*gr_jCoordsGlb(RIGHT_EDGE,gJ,1) + 4.0*gr_jCoordsGlb(RIGHT_EDGE,gJ-1,1) - gr_jCoordsGlb(RIGHT_EDGE,gJ-2,1))
    gr_jMetricsGlb(LEFT_EDGE, 2:gJ,1)   = 1.0 / (gr_jCoordsGlb(CENTER,2:gJ,1)     - gr_jCoordsGlb(CENTER,1:gJ-1,1))
    gr_jMetricsGlb(LEFT_EDGE, 1,1)      = 2.0 / (-3.0*gr_jCoordsGlb(LEFT_EDGE,1,1)   + 4.0*gr_jCoordsGlb(LEFT_EDGE,2,1)     - gr_jCoordsGlb(LEFT_EDGE,3,1))

  else ! 1d simulation

    gr_jCoords(LEFT_EDGE:RIGHT_EDGE,1,1)=gr_jmin
    gr_delta(JAXIS)=0.0

  end if

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Create z grid coordinates and metrics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(NDIM > 2) then ! 3d simulation

    if(gr_kStr) then !stretched grid in z
      select case (gr_kStrType)

      case(SG_USER)

        ! stretching method SG_USER not currently supported
        print*,"Create Domain : invalid stretching method "
        call Driver_abortFlash("Create Domain : invalid stretching method ")

      case(SG_TANH)

        ! start with a uniform distribution of points [0,1]
        gr_delta(KAXIS) = 1.0/gr_gIndexSize(KAXIS)
        halfDelta = gr_delta(KAXIS)/2.0
        j = gr_blkCornerID(KAXIS)-gr_guard(KAXIS)-1
        
        ! transform uniform to stretched using z(s) = Tanh((-1+2s) atan(a)) + 1 / 2a
        do i = gr_kloGc,gr_khiGc
          gr_kCoords(LEFT_EDGE,i,1)  = (gr_kmax-gr_kmin)*(tanh((-1.0+2.0*j*gr_delta(KAXIS)          )*atanh(gr_kStrPar))+1.0)/(2.0*gr_kStrPar)+gr_kmin
          gr_kCoords(CENTER,i,1)     = (gr_kmax-gr_kmin)*(tanh((-1.0+2.0*j*gr_delta(KAXIS)+halfDelta)*atanh(gr_kStrPar))+1.0)/(2.0*gr_kStrPar)+gr_kmin
          j = j+1
          gr_kCoords(RIGHT_EDGE,i,1) = (gr_kmax-gr_kmin)*(tanh((-1.0+2.0*j*gr_delta(KAXIS)          )*atanh(gr_kStrPar))+1.0)/(2.0*gr_kStrPar)+gr_kmin
        end do

        ! calculate global cell KAXIS coordinates
        j = 0
        do i = 1,gK
          gr_kCoordsGlb(LEFT_EDGE,i,1)  = (gr_kmax-gr_kmin)*(tanh((-1.0+2.0*j*gr_delta(KAXIS)          )*atanh(gr_kStrPar))+1.0)/(2.0*gr_kStrPar)+gr_kmin
          gr_kCoordsGlb(CENTER,i,1)     = (gr_kmax-gr_kmin)*(tanh((-1.0+2.0*j*gr_delta(KAXIS)+halfDelta)*atanh(gr_kStrPar))+1.0)/(2.0*gr_kStrPar)+gr_kmin
          j = j+1
          gr_kCoordsGlb(RIGHT_EDGE,i,1) = (gr_kmax-gr_kmin)*(tanh((-1.0+2.0*j*gr_delta(KAXIS)          )*atanh(gr_kStrPar))+1.0)/(2.0*gr_kStrPar)+gr_kmin
        end do
      
      case default

        ! stretching method unknown
        print*,"Create Domain : invalid stretching method "
        call Driver_abortFlash("Create Domain : invalid stretching method ")

      end select

    else ! uniform grid
	
      ! calculate block cell KAXIS coordinates      
      gr_delta(KAXIS) = (gr_kmax-gr_kmin)/gr_gIndexSize(KAXIS)
      halfDelta = gr_delta(KAXIS)/2.0
      j = gr_blkCornerID(KAXIS)-gr_guard(KAXIS)-1

      do i = gr_kloGc,gr_khiGc
        gr_kCoords(LEFT_EDGE,i,1) = gr_kmin+j*gr_delta(KAXIS)
        gr_kCoords(CENTER,i,1) = gr_kmin+j*gr_delta(KAXIS)+halfDelta
        j = j+1
        gr_kCoords(RIGHT_EDGE,i,1) = gr_kmin+j*gr_delta(KAXIS)
     end do

     ! calculate global cell KAXIS coordinates
     j = 0
     do i = 1,gK
       gr_kCoordsGlb(LEFT_EDGE,i,1) = gr_kmin+j*gr_delta(KAXIS)
       gr_kCoordsGlb(CENTER,i,1) = gr_kmin+j*gr_delta(KAXIS)+halfDelta
       j = j+1
       gr_kCoordsGlb(RIGHT_EDGE,i,1) = gr_kmin+j*gr_delta(KAXIS)
     end do

    end if    

    if(gr_axisMe(KAXIS)/=0)gr_blkBC(LOW,KAXIS)=NOT_BOUNDARY
    if(gr_axisMe(KAXIS)/=(gr_axisNumProcs(KAXIS)-1))&
          gr_blkBC(HIGH,KAXIS)=NOT_BOUNDARY

    ! calculate block KAXIS direction metric coefficients
    gr_kMetrics(CENTER,:,1)                       = 1.0 / (gr_kCoords(RIGHT_EDGE,gr_kloGC:gr_khiGC,1) - gr_kCoords(LEFT_EDGE,gr_kloGC:gr_khiGC,1))
    gr_kMetrics(RIGHT_EDGE,gr_kloGC:gr_khiGC-1,1) = 1.0 / (gr_kCoords(CENTER,gr_kloGC+1:gr_khiGC,1)   - gr_kCoords(CENTER,gr_kloGC:gr_khiGC-1,1))
    gr_kMetrics(RIGHT_EDGE,gr_khiGC,1)            = 2.0 / (-3.0*gr_kCoords(RIGHT_EDGE,gr_khiGC,1) + 4.0*gr_kCoords(RIGHT_EDGE,gr_khiGC-1,1) - gr_kCoords(RIGHT_EDGE,gr_khiGC-2,1))
    gr_kMetrics(LEFT_EDGE, gr_kloGC+1:gr_khiGC,1) = 1.0 / (gr_kCoords(CENTER,gr_kloGC+1:gr_khiGC,1)   - gr_kCoords(CENTER,gr_kloGC:gr_khiGC-1,1))
    gr_kMetrics(LEFT_EDGE, gr_kloGC,1)            = 2.0 / (-3.0*gr_kCoords(LEFT_EDGE,gr_kloGC,1)  + 4.0*gr_kCoords(LEFT_EDGE,gr_kloGC+1,1)  - gr_kCoords(LEFT_EDGE,gr_kloGC+2,1))

    ! calculate global KAXIS direction metric coefficients
    gr_kMetricsGlb(CENTER,1:gK,1)       = 1.0 / (gr_kCoordsGlb(RIGHT_EDGE,1:gK,1) - gr_kCoordsGlb(LEFT_EDGE,1:gK,1))
    gr_kMetricsGlb(RIGHT_EDGE,1:gK-1,1) = 1.0 / (gr_kCoordsGlb(CENTER,2:gK,1)     - gr_kCoordsGlb(CENTER,1:gK-1,1))
    gr_kMetricsGlb(RIGHT_EDGE,gK,1)     = 2.0 / (-3.0*gr_kCoordsGlb(RIGHT_EDGE,gK,1) + 4.0*gr_kCoordsGlb(RIGHT_EDGE,gK-1,1) - gr_kCoordsGlb(RIGHT_EDGE,gK-2,1))
    gr_kMetricsGlb(LEFT_EDGE, 2:gK,1)   = 1.0 / (gr_kCoordsGlb(CENTER,2:gK,1)     - gr_kCoordsGlb(CENTER,1:gK-1,1))
    gr_kMetricsGlb(LEFT_EDGE, 1,1)      = 2.0 / (-3.0*gr_kCoordsGlb(LEFT_EDGE,1,1)   + 4.0*gr_kCoordsGlb(LEFT_EDGE,2,1)     - gr_kCoordsGlb(LEFT_EDGE,3,1))

  else ! 2d simulation

    gr_kCoords(LEFT_EDGE:RIGHT_EDGE,1,1) = gr_kmin
    gr_delta(KAXIS)=0.0

  end if
 
  deallocate(gr_iCoordsGlb)
  deallocate(gr_jCoordsGlb)
  deallocate(gr_kCoordsGlb) 
  
end subroutine gr_createDomain
