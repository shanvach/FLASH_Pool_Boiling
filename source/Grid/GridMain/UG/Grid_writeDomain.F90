!!****if* source/Grid/GridMain/UG/Grid_writeDomain
!!
!! NAME
!!
!!  Grid_writeDomain
!!
!!
!! SYNOPSIS
!!
!!  call Grid_writeDomain()
!!
!!
!! DESCRIPTION
!!
!! This function writes the grid information to an hdf5 file to store the 
!! Uniform Grid / Regular Grid cell coordinates (Left, Center, Right) and 
!! the cell metrics for later use in post-processing FLASH simulations.
!!
!! Currently only supports hdf5 IO
!!
!! ARGUMENTS
!!
!!  none
!!
!!***
 
subroutine Grid_writeDomain()

  use Grid_data, ONLY : gr_meshMe, gr_meshComm, gr_globalNumBlocks,    &
                        gr_iCoords, gr_jCoords, gr_kCoords,            &
                        gr_iMetrics, gr_jMetrics, gr_kMetrics,         &
                        gr_ilo, gr_ihi, gr_jlo, gr_jhi, gr_klo, gr_khi
  use Timers_interface, ONLY : Timers_start, Timers_stop

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"

  integer :: fileID, ierr, i, j, k, a, b, c, &
             realSize, newType
  integer(kind=MPI_ADDRESS_KIND) :: extent, begin
  integer, dimension(2) :: sizes, subSizes, starts
  integer, dimension(MDIM) :: axes, resizedType
  integer, allocatable, dimension(:) :: counts, displs
  character(len=MAX_STRING_LENGTH) :: filename
  character(len=5) :: gCrdLbs(MDIM, 3), gMtrLbs(MDIM, 3)
  real, dimension(MDIM, 3) :: gCrdMax, gCrdMin, gMtrMax, gMtrMin 
  real, allocatable, dimension(:,:) :: iCoords, jCoords, kCoords
  real, allocatable, dimension(:,:) :: iMetrics, jMetrics, kMetrics
  real, allocatable, dimension(:,:,:,:,:,:,:) :: gCoords, gMetrics

#ifdef FLASH_IO_HDF5

  call Timers_start("writeDomain")

  ! Open an hdf5 file for writing grid information
  if(gr_meshMe == MASTER_PE) then

    call io_getOutputName(0, "hdf5_", "grd_", filename, .false.)

    fileID = -1
    call io_h5init_file(fileID, filename)
    if(fileID == -1) then
      print *, "Error: Unable to initialize GRID OUTPUT file"
      call Driver_abortFlash("Unable to initialize hdf5 file")
    end if

  end if
  
  ! Create mpi subdomain parameters (const across axes)
  allocate(counts(gr_globalNumBlocks), displs(gr_globalNumBlocks))
  begin = 0
  counts = 1
  starts = [0, 0]
  call MPI_Type_size(FLASH_REAL, realSize, ierr)
  extent = 3 * realSize
  do k=1, gr_globalNumBlocks
    displs(k) = (k-1) 
  end do

  ! Create MPI TYPE for comms for each axis
  axes = (/ NXB, NYB, NZB /)
  do k=IAXIS, KAXIS
    sizes = [3*gr_globalNumBlocks, axes(k)] 
    subSizes = [3, axes(k)] 
    call MPI_Type_create_subarray(2, sizes, subSizes, starts,    &
                                  MPI_ORDER_FORTRAN, FLASH_REAL, &
                                  newType, ierr)
    call MPI_Type_create_resized(newType, begin, extent, resizedType(k), ierr)
    call MPI_Type_commit(resizedType(k), ierr)
  end do

  if(gr_meshMe == MASTER_PE) then

    ! Create dataset labels
    gCrdLbs(IAXIS, LEFT_EDGE:RIGHT_EDGE) = (/ 'xxxl', 'xxxc', 'xxxr' /) 
    gCrdLbs(JAXIS, LEFT_EDGE:RIGHT_EDGE) = (/ 'yyyl', 'yyyc', 'yyyr' /) 
    gCrdLbs(KAXIS, LEFT_EDGE:RIGHT_EDGE) = (/ 'zzzl', 'zzzc', 'zzzr' /) 
    gMtrLbs(IAXIS, LEFT_EDGE:RIGHT_EDGE) = (/ 'ddxl', 'ddxc', 'ddxr' /) 
    gMtrLbs(JAXIS, LEFT_EDGE:RIGHT_EDGE) = (/ 'ddyl', 'ddyc', 'ddyr' /) 
    gMtrLbs(KAXIS, LEFT_EDGE:RIGHT_EDGE) = (/ 'ddzl', 'ddzc', 'ddzr' /) 

    ! Create mpi buffers and global grid storage array
    !   buffers shape are    (faces * blks, blk size)
    !   global grid shape is (axes, faces, blks, bsz I, bsz J, bsz K, 1)
    allocate(iCoords(3 * gr_globalNumBlocks, NXB))
    allocate(jCoords(3 * gr_globalNumBlocks, NYB))
    allocate(kCoords(3 * gr_globalNumBlocks, NZB))
    allocate(iMetrics(3 * gr_globalNumBlocks, NXB))
    allocate(jMetrics(3 * gr_globalNumBlocks, NYB))
    allocate(kMetrics(3 * gr_globalNumBlocks, NZB))
    allocate(gCoords(MDIM, 3, gr_globalNumBlocks, NXB, NYB, NZB, 1))
    allocate(gMetrics(MDIM, 3, gr_globalNumBlocks, NXB, NYB, NZB, 1))

  endif  
  
  ! Gather grid information from mpi from all processors
  call MPI_Gatherv(gr_iCoords(:,gr_ilo:gr_ihi,1), 3*NXB, FLASH_REAL, iCoords, &
                   counts, displs, resizedType(IAXIS), 0, gr_meshComm, ierr)
  call MPI_Gatherv(gr_jCoords(:,gr_jlo:gr_jhi,1), 3*NYB, FLASH_REAL, jCoords, &
                   counts, displs, resizedType(JAXIS), 0, gr_meshComm, ierr)
  call MPI_Gatherv(gr_kCoords(:,gr_klo:gr_khi,1), 3*NZB, FLASH_REAL, kCoords, &
                   counts, displs, resizedType(KAXIS), 0, gr_meshComm, ierr)
  call MPI_Gatherv(gr_iMetrics(:,gr_ilo:gr_ihi,1), 3*NXB, FLASH_REAL, iMetrics, &
                   counts, displs, resizedType(IAXIS), 0, gr_meshComm, ierr)
  call MPI_Gatherv(gr_jMetrics(:,gr_jlo:gr_jhi,1), 3*NYB, FLASH_REAL, jMetrics, &
                   counts, displs, resizedType(JAXIS), 0, gr_meshComm, ierr)
  call MPI_Gatherv(gr_kMetrics(:,gr_klo:gr_khi,1), 3*NZB, FLASH_REAL, kMetrics, &
                   counts, displs, resizedType(KAXIS), 0, gr_meshComm, ierr)
  deallocate(counts, displs)

  if(gr_meshMe == MASTER_PE) then

    ! Copy mpi buffers to global grid storage array
    gCoords(IAXIS,:,:,:,1,1,1) = reshape(iCoords, (/ 3, gr_globalNumBlocks, NXB /)) 
    gCoords(JAXIS,:,:,1,:,1,1) = reshape(jCoords, (/ 3, gr_globalNumBlocks, NYB /)) 
    gCoords(KAXIS,:,:,1,1,:,1) = reshape(kCoords, (/ 3, gr_globalNumBlocks, NZB /)) 
    gMetrics(IAXIS,:,:,:,1,1,1) = reshape(iMetrics, (/ 3, gr_globalNumBlocks, NXB /)) 
    gMetrics(JAXIS,:,:,1,:,1,1) = reshape(jMetrics, (/ 3, gr_globalNumBlocks, NYB /)) 
    gMetrics(KAXIS,:,:,1,1,:,1) = reshape(kMetrics, (/ 3, gr_globalNumBlocks, NZB /)) 
    deallocate(iCoords, jCoords, kCoords, iMetrics, jMetrics, kMetrics)

    ! fill out mesh grid
#if NDIM == 3
      do k=2, NZB
        do j=2, NYB
          gCoords(IAXIS,:,:,:,j,k,1) = gCoords(IAXIS,CENTER,:,:,j-1,k-1,1)
          gMetrics(IAXIS,:,:,:,j,k,1) = gMetrics(IAXIS,CENTER,:,:,j-1,k-1,1)
        end do
        do i=2, NXB
          gCoords(JAXIS,:,:,i,:,k,1) = gCoords(JAXIS,:,:,i-1,:,k-1,1)
          gMetrics(JAXIS,:,:,i,:,k,1) = gMetrics(JAXIS,:,:,i-1,:,k-1,1)
        end do
      end do
      do j=2, NYB
        do i=2, NXB
          gCoords(KAXIS,:,:,i,j,:,1) = gCoords(KAXIS,:,:,i-1,j-1,:,1)
          gMetrics(KAXIS,:,:,i,j,:,1) = gMetrics(KAXIS,:,:,i-1,j-1,:,1)
        end do
      end do
#else
      do j=2, NYB
        gCoords(IAXIS,:,:,:,j,1,1) = gCoords(IAXIS,:,:,:,j-1,1,1)
        gMetrics(IAXIS,:,:,:,j,1,1) = gMetrics(IAXIS,:,:,:,j-1,1,1)
      end do
      do i=2, NXB
        gCoords(JAXIS,:,:,i,:,1,1) = gCoords(JAXIS,:,:,i-1,:,1,1)
        gMetrics(JAXIS,:,:,i,:,1,1) = gMetrics(JAXIS,:,:,i-1,:,1,1)
      end do
      gCoords(KAXIS,:,:,:,:,1,1) = gCoords(KAXIS,1,1,1,1,1,1) 
      gMetrics(KAXIS,:,:,:,:,1,1) = gMetrics(KAXIS,1,1,1,1,1,1)
#endif      

    ! Write IAXIS coordinates 
    do a=IAXIS, KAXIS
      do b=LEFT_EDGE, RIGHT_EDGE
        do c=1, gr_globalNumBlocks    

          ! find extreme values
          gCrdMax(a, b) = maxval(gCoords(a, b, :, :, :, :, :)) 
          gCrdMax(a, b) = maxval(gCoords(a, b, :, :, :, :, :)) 
          gMtrMax(a, b) = maxval(gMetrics(a, b, :, :, :, :, :)) 
          gMtrMax(a, b) = maxval(gMetrics(a, b, :, :, :, :, :)) 
      
          ! write to hdf5 file
          call io_h5write_unknowns(gr_meshMe, fileID, NXB, NYB, NZB, & 
                                   gCoords(a, b, c, :, :, :, :),     &
                                   gCrdMin(a, b),                    & 
                                   gCrdMax(a, b),                    &
                                   gCrdLbs(a, b),                    &
                                   1, gr_globalNumBlocks, c-1)   
          call io_h5write_unknowns(gr_meshMe, fileID, NXB, NYB, NZB, & 
                                   gMetrics(a, b, c, :, :, :, :),    &
                                   gMtrMin(a, b),                    & 
                                   gMtrMax(a, b),                    &
                                   gMtrLbs(a, b),                    &
                                   1, gr_globalNumBlocks, c-1)   
        end do
      end do
    end do
   
    ! release storage arrays
    deallocate(gCoords, gMetrics)

  endif

  ! free MPI Type for comms
  call MPI_Type_Free(newType, ierr)
  do k=IAXIS, KAXIS
    call MPI_Type_Free(resizedType(k), ierr)
  end do

  ! Close the hdf5 file
  if(gr_meshMe == MASTER_PE) then

    call io_h5close_file(fileID)

  end if

  call Timers_stop("writeDomain")

#endif

  return

end subroutine Grid_WriteDomain
