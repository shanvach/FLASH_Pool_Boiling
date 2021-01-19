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
!! Paramesh or Uniform Grid / Regular Grid cell coordinates (Left, Center, Right)
!! and the cell metrics for later use in post-processing FLASH simulations.
!!
!! Currently only supports hdf5 IO
!!
!! This fuction is intended to be called after a IO_output function
!!
!! ARGUMENTS
!!
!!  none
!!
!!***
 
subroutine Grid_writeDomain(fileNumber)

  use Driver_interface, ONLY : Driver_abortFlash
  use IO_data, ONLY : io_plotFileNumber
  use Grid_data, ONLY : gr_meshMe, gr_meshComm, gr_globalNumBlocks,    &
                        gr_iCoords, gr_jCoords, gr_kCoords, gr_delta,  &
                        gr_ilo, gr_ihi, gr_jlo, gr_jhi, gr_klo, gr_khi
#ifdef FLASH_GRID_REGULAR
  use Grid_data, ONLY : gr_iMetrics, gr_jMetrics, gr_kMetrics
#endif

  use Timers_interface, ONLY : Timers_start, Timers_stop

  use HDF5

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"

  integer, optional, intent(IN) :: fileNumber

  integer :: ierr, i, j, k, a, b, c, d, &
             realSize, newType, jproc
  integer(kind=MPI_ADDRESS_KIND) :: extent, begin
  integer :: status(MPI_STATUS_SIZE)
  integer, dimension(2) :: sizes, subSizes, starts
  integer, dimension(MDIM) :: axes, resizedType
  integer, allocatable, dimension(:) :: counts, displs
  character(len=5) :: gCrdLbs(MDIM, 4), gMtrLbs(MDIM, 3)
  logical :: force_gatherv, used_gatherv
  real, dimension(MDIM, 4) :: gCrdMax, gCrdMin, gMtrMax, gMtrMin 
  real, allocatable, dimension(:,:) :: iCoords, jCoords, kCoords
  real, allocatable, dimension(:) :: iBuff, jBuff, kBuff 
  real, allocatable, dimension(:,:) :: fCoord
  real, allocatable, dimension(:,:) :: iMetrics, jMetrics, kMetrics
  real, allocatable, dimension(:,:,:,:,:) :: gCoords, gMetrics
  real, allocatable, dimension(:,:) :: gDeltas

  ! locals necessary to read hdf5 file
  integer :: error
  integer(HID_T) :: file_id, dspc_id, dset_id, aspc_id, attr_id 
  integer(HSIZE_T), dimension(2) :: dset_dims
  integer(HSIZE_T), dimension(1) :: dset_sngl
  character(len = 32) :: dsetname, attrname
  character(len = MAX_STRING_LENGTH) :: filename

#ifdef FLASH_IO_HDF5

  call Timers_start("writeDomain")

  ! Open an hdf5 file for writing grid information
  if(gr_meshMe == MASTER_PE) then

    if (present(filenumber)) then
      call io_getOutputName(fileNumber, "hdf5_", "grd_", filename, .false.)
    else
      call io_getOutputName(max(io_plotFileNumber - 1, 0), "hdf5_", "grd_", filename, .false.)
    endif

    call h5open_f(error)

    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
    if(file_id == -1) then
      print *, "Error: Unable to initialize GRID OUTPUT file"
      call Driver_abortFlash("Unable to initialize hdf5 file")
    end if

  end if
  
  if(gr_meshMe == MASTER_PE) then

    ! Create dataset labels
    gCrdLbs(IAXIS, LEFT_EDGE:RIGHT_EDGE + 1) = (/ 'xxxl', 'xxxc', 'xxxr', 'xxxf' /) 
    gCrdLbs(JAXIS, LEFT_EDGE:RIGHT_EDGE + 1) = (/ 'yyyl', 'yyyc', 'yyyr', 'yyyf' /) 
    gCrdLbs(KAXIS, LEFT_EDGE:RIGHT_EDGE + 1) = (/ 'zzzl', 'zzzc', 'zzzr', 'zzzf' /) 
    gMtrLbs(IAXIS, LEFT_EDGE:RIGHT_EDGE) = (/ 'ddxl', 'ddxc', 'ddxr' /) 
    gMtrLbs(JAXIS, LEFT_EDGE:RIGHT_EDGE) = (/ 'ddyl', 'ddyc', 'ddyr' /) 
    gMtrLbs(KAXIS, LEFT_EDGE:RIGHT_EDGE) = (/ 'ddzl', 'ddzc', 'ddzr' /) 

    ! Create mpi buffers and global grid storage array
    !   buffers shape are    (faces * blks, blk size)
    !   global grid shape is (axes, faces, blks, bsz, 1)
    allocate(iCoords(3 * gr_globalNumBlocks, NXB))
    allocate(jCoords(3 * gr_globalNumBlocks, NYB))
    allocate(kCoords(3 * gr_globalNumBlocks, NZB))
    allocate(gCoords(MDIM, 3, gr_globalNumBlocks, max(NXB, NYB, NZB), 1))

#ifdef FLASH_GRID_REGULAR
    allocate(iMetrics(3 * gr_globalNumBlocks, NXB))
    allocate(jMetrics(3 * gr_globalNumBlocks, NYB))
    allocate(kMetrics(3 * gr_globalNumBlocks, NZB))
    allocate(gMetrics(MDIM, 3, gr_globalNumBlocks, max(NXB, NYB, NZB), 1))
#else
    allocate(gDeltas(gr_globalNumBlocks, MDIM))
#endif

  endif  

  ! Sometimes there is an issue if the gatherv buffer is too large
  force_gatherv = .false.
  if (NXB <= 256 .and. NYB <= 256 .and. NZB <= 256 .and. gr_globalNumBlocks <= 64 .or. force_gatherv) then
    used_gatherv = .true.

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
    
    ! Gather grid information from mpi from all processors
    call MPI_Gatherv(gr_iCoords(:,gr_ilo:gr_ihi,1), 3*NXB, FLASH_REAL, iCoords, &
                     counts, displs, resizedType(IAXIS), 0, gr_meshComm, ierr)
    call MPI_Gatherv(gr_jCoords(:,gr_jlo:gr_jhi,1), 3*NYB, FLASH_REAL, jCoords, &
                     counts, displs, resizedType(JAXIS), 0, gr_meshComm, ierr)
    call MPI_Gatherv(gr_kCoords(:,gr_klo:gr_khi,1), 3*NZB, FLASH_REAL, kCoords, &
                     counts, displs, resizedType(KAXIS), 0, gr_meshComm, ierr)

#ifdef FLASH_GRID_REGULAR
    call MPI_Gatherv(gr_iMetrics(:,gr_ilo:gr_ihi,1), 3*NXB, FLASH_REAL, iMetrics, &
                     counts, displs, resizedType(IAXIS), 0, gr_meshComm, ierr)
    call MPI_Gatherv(gr_jMetrics(:,gr_jlo:gr_jhi,1), 3*NYB, FLASH_REAL, jMetrics, &
                     counts, displs, resizedType(JAXIS), 0, gr_meshComm, ierr)
    call MPI_Gatherv(gr_kMetrics(:,gr_klo:gr_khi,1), 3*NZB, FLASH_REAL, kMetrics, &
                     counts, displs, resizedType(KAXIS), 0, gr_meshComm, ierr)
#else
    call MPI_Gather(gr_delta(IAXIS), 1, FLASH_REAL, gDeltas(:, IAXIS), &
                                     1, FLASH_REAL, 0, gr_meshComm, ierr)
    call MPI_Gather(gr_delta(JAXIS), 1, FLASH_REAL, gDeltas(:, JAXIS), &
                                     1, FLASH_REAL, 0, gr_meshComm, ierr)
    call MPI_Gather(gr_delta(KAXIS), 1, FLASH_REAL, gDeltas(:, KAXIS), &
                                     1, FLASH_REAL, 0, gr_meshComm, ierr)
#endif

    deallocate(counts, displs)

  ! fall back to point to point communications
  else
    used_gatherv = .false.

    do jproc = 0, gr_globalNumBlocks - 1

      allocate(iBuff(3*NXB), jBuff(3*NYB), kBuff(3*NZB))

      if (gr_meshMe == MASTER_PE .and. jproc == MASTER_PE) then
        iCoords(3*jproc+1:3*jproc+3,:) = gr_iCoords(:,gr_ilo:gr_ihi,1)
        jCoords(3*jproc+1:3*jproc+3,:) = gr_jCoords(:,gr_jlo:gr_jhi,1)
        kCoords(3*jproc+1:3*jproc+3,:) = gr_kCoords(:,gr_klo:gr_khi,1)

#ifdef FLASH_GRID_REGULAR
        iMetrics(3*jproc+1:3*jproc+3,:) = gr_iMetrics(:,gr_ilo:gr_ihi,1)
        jMetrics(3*jproc+1:3*jproc+3,:) = gr_jMetrics(:,gr_jlo:gr_jhi,1)
        kMetrics(3*jproc+1:3*jproc+3,:) = gr_kMetrics(:,gr_klo:gr_khi,1)
#else
        gDeltas(jproc, IAXIS) = gr_delta(IAXIS)
        gDeltas(jproc, JAXIS) = gr_delta(JAXIS)
        gDeltas(jproc, KAXIS) = gr_delta(KAXIS)
#endif
      endif

      if (gr_meshMe == MASTER_PE .and. jproc /= MASTER_PE) then
        call MPI_Recv(iBuff, 3*NXB, FLASH_REAL, jproc, IAXIS+1, gr_meshComm, status, ierr)
        call MPI_Recv(jBuff, 3*NYB, FLASH_REAL, jproc, JAXIS+1, gr_meshComm, status, ierr)
        call MPI_Recv(kBuff, 3*NZB, FLASH_REAL, jproc, KAXIS+1, gr_meshComm, status, ierr)
       
        iCoords(3*jproc+1:3*jproc+3,:) = reshape(iBuff, (/ 3, NXB /))
        jCoords(3*jproc+1:3*jproc+3,:) = reshape(jBuff, (/ 3, NYB /))
        kCoords(3*jproc+1:3*jproc+3,:) = reshape(kBuff, (/ 3, NZB /))

#ifdef FLASH_GRID_REGULAR
        call MPI_Recv(iBuff, 3*NXB, FLASH_REAL, jproc, IAXIS+2, gr_meshComm, status, ierr)
        call MPI_Recv(jBuff, 3*NYB, FLASH_REAL, jproc, JAXIS+2, gr_meshComm, status, ierr)
        call MPI_Recv(kBuff, 3*NZB, FLASH_REAL, jproc, KAXIS+2, gr_meshComm, status, ierr)
       
        iMetrics(3*jproc+1:3*jproc+3,:) = reshape(iBuff, (/ 3, NXB /))
        jMetrics(3*jproc+1:3*jproc+3,:) = reshape(jBuff, (/ 3, NYB /))
        kMetrics(3*jproc+1:3*jproc+3,:) = reshape(kBuff, (/ 3, NZB /))
#else
        call MPI_Recv(gDeltas(jproc, IAXIS), 1, FLASH_REAL, jproc, IAXIS+2, gr_meshComm, status, ierr)
        call MPI_Recv(gDeltas(jproc, JAXIS), 1, FLASH_REAL, jproc, JAXIS+2, gr_meshComm, status, ierr)
        call MPI_Recv(gDeltas(jproc, KAXIS), 1, FLASH_REAL, jproc, KAXIS+2, gr_meshComm, status, ierr)
#endif
      endif

      if (gr_meshMe /= MASTER_PE .and. jproc == gr_meshMe) then
        iBuff = reshape(gr_iCoords(1:3,gr_ilo:gr_ihi,1), (/ 3*NXB /))
        jBuff = reshape(gr_jCoords(1:3,gr_jlo:gr_jhi,1), (/ 3*NYB /))
        kBuff = reshape(gr_kCoords(1:3,gr_klo:gr_khi,1), (/ 3*NZB /))

        call MPI_Send(iBuff, 3*NXB, FLASH_REAL, MASTER_PE, IAXIS+1, gr_meshComm, ierr)
        call MPI_Send(jBuff, 3*NYB, FLASH_REAL, MASTER_PE, JAXIS+1, gr_meshComm, ierr)
        call MPI_Send(kBuff, 3*NZB, FLASH_REAL, MASTER_PE, KAXIS+1, gr_meshComm, ierr)

#ifdef FLASH_GRID_REGULAR
        iBuff = reshape(gr_iMetrics(1:3,gr_ilo:gr_ihi,1), (/ 3*NXB /))
        jBuff = reshape(gr_jMetrics(1:3,gr_jlo:gr_jhi,1), (/ 3*NYB /))
        kBuff = reshape(gr_kMetrics(1:3,gr_klo:gr_khi,1), (/ 3*NZB /))

        call MPI_Send(iBuff, 3*NXB, FLASH_REAL, MASTER_PE, IAXIS+2, gr_meshComm, ierr)
        call MPI_Send(jBuff, 3*NYB, FLASH_REAL, MASTER_PE, JAXIS+2, gr_meshComm, ierr)
        call MPI_Send(kBuff, 3*NZB, FLASH_REAL, MASTER_PE, KAXIS+2, gr_meshComm, ierr)
#else
        call MPI_Send(gr_delta(IAXIS), 1, FLASH_REAL, MASTER_PE, IAXIS+2, gr_meshComm, ierr)
        call MPI_Send(gr_delta(JAXIS), 1, FLASH_REAL, MASTER_PE, JAXIS+2, gr_meshComm, ierr)
        call MPI_Send(gr_delta(KAXIS), 1, FLASH_REAL, MASTER_PE, KAXIS+2, gr_meshComm, ierr)
#endif
      endif

      deallocate(iBuff, jBuff, kBuff)

    end do

  endif

  if(gr_meshMe == MASTER_PE) then

    ! Copy mpi buffers to global grid storage array
    gCoords(IAXIS,:,:,1:NXB,1) = reshape(iCoords, (/ 3, gr_globalNumBlocks, NXB /)) 
    gCoords(JAXIS,:,:,1:NYB,1) = reshape(jCoords, (/ 3, gr_globalNumBlocks, NYB /)) 
    gCoords(KAXIS,:,:,1:NZB,1) = reshape(kCoords, (/ 3, gr_globalNumBlocks, NZB /)) 
    deallocate(iCoords, jCoords, kCoords)
#ifdef FLASH_GRID_REGULAR
    gMetrics(IAXIS,:,:,1:NXB,1) = reshape(iMetrics, (/ 3, gr_globalNumBlocks, NXB /)) 
    gMetrics(JAXIS,:,:,1:NYB,1) = reshape(jMetrics, (/ 3, gr_globalNumBlocks, NYB /)) 
    gMetrics(KAXIS,:,:,1:NZB,1) = reshape(kMetrics, (/ 3, gr_globalNumBlocks, NZB /)) 
    deallocate(iMetrics, jMetrics, kMetrics)
#endif

    ! Write coordinates to file
    do a=IAXIS, KAXIS 

      ! Determine bounds      
      select case (a)
        case (IAXIS)
          d = NXB
        case (JAXIS)
          d = NYB
        case (KAXIS)
          d = NZB
      end select 

      ! write each face per axis
      do b=LEFT_EDGE, RIGHT_EDGE

        ! find extreme values
        gCrdMax(a, b) = maxval(gCoords(a, b, :, 1:d, :)) 
        gCrdMin(a, b) = minval(gCoords(a, b, :, 1:d, :)) 
#ifdef FLASH_GRID_REGULAR
        gMtrMax(a, b) = maxval(gMetrics(a, b, :, 1:d, :)) 
        gMtrMin(a, b) = minval(gMetrics(a, b, :, 1:d, :)) 
#else
        gMtrMax(a, b) = maxval(gDeltas(:, a)) 
        gMtrMin(a, b) = minval(gDeltas(:, a)) 
#endif

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
#ifdef FLASH_GRID_REGULAR
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
#else
        if (b == CENTER) then
          dsetname = gMtrLbs(a, b)
          dset_sngl = (/ gr_globalNumBlocks/)
          call h5screate_simple_f(1, dset_sngl, dspc_id, error)
          call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspc_id, dset_id, error)  
          call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, gDeltas(:, a), dset_dims, error)
       
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
        endif
#endif

      end do

      ! Consolodate axis "face" points
      allocate(fCoord(gr_globalNumBlocks, d + 1))
      fCoord(:, 1) = gCoords(a, LEFT_EDGE, :, 1, 1)
      fCoord(:, 2:d+1) = gCoords(a, RIGHT_EDGE, :, 1:d, 1) 
      if (a == KAXIS .AND. d == 1) fCoord(:, 2) = fCoord(:, 1) + 0.000001

      ! find extreme values
      gCrdMax(a, 4) = maxval(fCoord) 
      gCrdMin(a, 4) = minval(fCoord) 
   
      ! write dimensions
      dsetname = gCrdLbs(a, 4)
      dset_dims = (/ d + 1, gr_globalNumBlocks /)
      call h5screate_simple_f(2, dset_dims, dspc_id, error)
      call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspc_id, dset_id, error)  
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, transpose(fCoord), dset_dims, error)
       
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

      deallocate(fCoord)

    end do
   
    ! release storage arrays
    deallocate(gCoords)
#ifdef FLASH_GRID_REGULAR
    deallocate(gMetrics)
#else
    deallocate(gDeltas)
#endif
  endif

  ! clean up collective comms
  if (used_gatherv) then
    call MPI_Type_Free(newType, ierr)
    do k=IAXIS, KAXIS
      call MPI_Type_Free(resizedType(k), ierr)
    end do
  endif

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
