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

  use Driver_interface, ONLY : Driver_abortFlash

  use Grid_data, ONLY : gr_meshMe, gr_meshComm, gr_globalNumBlocks,    &
                        gr_iCoords, gr_jCoords, gr_kCoords,            &
                        gr_iMetrics, gr_jMetrics, gr_kMetrics,         &
                        gr_ilo, gr_ihi, gr_jlo, gr_jhi, gr_klo, gr_khi

  use Timers_interface, ONLY : Timers_start, Timers_stop

  use HDF5

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"

  integer :: ierr, i, j, k, a, b, c, d, &
             realSize, newType
  integer(kind=MPI_ADDRESS_KIND) :: extent, begin
  integer, dimension(2) :: sizes, subSizes, starts
  integer, dimension(MDIM) :: axes, resizedType
  integer, allocatable, dimension(:) :: counts, displs
  character(len=5) :: gCrdLbs(MDIM, 3), gMtrLbs(MDIM, 3)
  real, dimension(MDIM, 3) :: gCrdMax, gCrdMin, gMtrMax, gMtrMin 
  real, allocatable, dimension(:,:) :: iCoords, jCoords, kCoords
  real, allocatable, dimension(:,:) :: iMetrics, jMetrics, kMetrics
  real, allocatable, dimension(:,:,:,:,:) :: gCoords, gMetrics

  ! locals necessary to read hdf5 file
  integer :: error
  integer(HID_T) :: file_id, dspc_id, dset_id, aspc_id, attr_id 
  integer(HSIZE_T), dimension(2) :: dset_dims
  integer(HSIZE_T), dimension(1) :: dset_sngl
  character(len = 32) :: filename, dsetname, attrname

#ifdef FLASH_IO_HDF5

  call Timers_start("writeDomain")


  ! Open an hdf5 file for writing grid information
  if(gr_meshMe == MASTER_PE) then

    call io_getOutputName(0, "hdf5_", "grd_", filename, .false.)

    call h5open_f(error)

    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
    if(file_id == -1) then
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
    !   global grid shape is (axes, faces, blks, bsz, 1)
    allocate(iCoords(3 * gr_globalNumBlocks, NXB))
    allocate(jCoords(3 * gr_globalNumBlocks, NYB))
    allocate(kCoords(3 * gr_globalNumBlocks, NZB))
    allocate(iMetrics(3 * gr_globalNumBlocks, NXB))
    allocate(jMetrics(3 * gr_globalNumBlocks, NYB))
    allocate(kMetrics(3 * gr_globalNumBlocks, NZB))
    allocate(gCoords(MDIM, 3, gr_globalNumBlocks, max(NXB, NYB, NZB), 1))
    allocate(gMetrics(MDIM, 3, gr_globalNumBlocks, max(NXB, NYB, NZB), 1))

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
    gCoords(IAXIS,:,:,1:NXB,1) = reshape(iCoords, (/ 3, gr_globalNumBlocks, NXB /)) 
    gCoords(JAXIS,:,:,1:NYB,1) = reshape(jCoords, (/ 3, gr_globalNumBlocks, NYB /)) 
    gCoords(KAXIS,:,:,1:NZB,1) = reshape(kCoords, (/ 3, gr_globalNumBlocks, NZB /)) 
    gMetrics(IAXIS,:,:,1:NXB,1) = reshape(iMetrics, (/ 3, gr_globalNumBlocks, NXB /)) 
    gMetrics(JAXIS,:,:,1:NYB,1) = reshape(jMetrics, (/ 3, gr_globalNumBlocks, NYB /)) 
    gMetrics(KAXIS,:,:,1:NZB,1) = reshape(kMetrics, (/ 3, gr_globalNumBlocks, NZB /)) 
    deallocate(iCoords, jCoords, kCoords, iMetrics, jMetrics, kMetrics)

    ! Write coordinates to file
    do a=IAXIS, KAXIS 

      ! Determine bounds      
      select case (a)
        case (1)
          d = NXB
        case (2)
          d = NYB
        case (3)
          d = NZB
      end select 

      ! write each face per axis
      do b=LEFT_EDGE, RIGHT_EDGE

        ! find extreme values
        gCrdMax(a, b) = maxval(gCoords(a, b, :, 1:d, :)) 
        gCrdMax(a, b) = maxval(gCoords(a, b, :, 1:d, :)) 
        gMtrMax(a, b) = maxval(gMetrics(a, b, :, 1:d, :)) 
        gMtrMax(a, b) = maxval(gMetrics(a, b, :, 1:d, :)) 
   
        ! write dimensions
        dsetname = gCrdLbs(a, b)
        dset_dims = (/ d, gr_globalNumBlocks /)
        call h5screate_simple_f(2, dset_dims, dspc_id, error)
        call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspc_id, dset_id, error)  
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, transpose(gCoords(a, b, :, 1:d, 1)), dset_dims, error)
       
        attrname = "maximum"
        dset_sngl = (/ 1 /)
        call h5screate_simple_f(1, dset_sngl, aspc_id, error)
        call h5acreate_f(dset_id, attrname, H5T_NATIVE_DOUBLE, aspc_id, attr_id, error) 
        call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, gCrdMax(a, b), dset_sngl, error)
        call h5aclose_f(attr_id, error)

        attrname = "minimum"
        call h5acreate_f(dset_id, attrname, H5T_NATIVE_DOUBLE, aspc_id, attr_id, error) 
        call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, gCrdMin(a, b), dset_sngl, error)
        call h5aclose_f(attr_id, error)
        call h5sclose_f(aspc_id, error)

        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspc_id, error)

        ! write metrics 
        dsetname = gMtrLbs(a, b)
        dset_dims = (/ d, gr_globalNumBlocks /)
        call h5screate_simple_f(2, dset_dims, dspc_id, error)
        call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspc_id, dset_id, error)  
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, transpose(gMetrics(a, b, :, 1:d, 1)), dset_dims, error)
       
        attrname = "maximum"
        dset_sngl = (/ 1 /)
        call h5screate_simple_f(1, dset_sngl, aspc_id, error)
        call h5acreate_f(dset_id, attrname, H5T_NATIVE_DOUBLE, aspc_id, attr_id, error) 
        call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, gMtrMax(a, b), dset_sngl, error)
        call h5aclose_f(attr_id, error)

        attrname = "minimum"
        call h5acreate_f(dset_id, attrname, H5T_NATIVE_DOUBLE, aspc_id, attr_id, error) 
        call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, gMtrMin(a, b), dset_sngl, error)
        call h5aclose_f(attr_id, error)
        call h5sclose_f(aspc_id, error)

        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspc_id, error)

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

    print *, "*** Wrote grid data to ", trim(filename), " ****"

    call h5fclose_f(file_id, error)
    call h5close_f(error)

  end if

  call Timers_stop("writeDomain")

#endif

  return

end subroutine Grid_WriteDomain
