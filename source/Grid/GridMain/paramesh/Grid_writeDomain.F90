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
  use tree, ONLY : lrefine, lnblocks
  use Grid_data, ONLY : gr_meshMe, gr_meshComm, gr_globalNumBlocks,    &
                        gr_oneBlock, gr_delta, gr_globalNumProcs,      &
                        gr_ilo, gr_ihi, gr_jlo, gr_jhi, gr_klo, gr_khi

  use Timers_interface, ONLY : Timers_start, Timers_stop

  use HDF5

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"

  integer :: ierr, fc, lb, a, b, c, d, jproc, offset, localNumBlocks
  integer :: status(MPI_STATUS_SIZE)
  character(len=5) :: gCrdLbs(MDIM, 4), gMtrLbs(MDIM, 3)
  real, dimension(MDIM, 4) :: gCrdMax, gCrdMin, gMtrMax, gMtrMin 
  real, allocatable, dimension(:) :: iBuff, jBuff, kBuff, dBuff 
  real, allocatable, dimension(:,:) :: fCoord
  real, allocatable, dimension(:,:,:,:) :: gCoords
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

    call io_getOutputName(0, "hdf5_", "grd_", filename, .false.)

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

    allocate(gCoords(max(NXB, NYB, NZB), gr_globalNumBlocks, 3, MDIM))
    allocate(gDeltas(gr_globalNumBlocks, MDIM))

  endif  

  offset = 0
  do jproc = 0, gr_globalNumProcs - 1


    if (gr_meshMe == MASTER_PE .and. jproc == MASTER_PE) then
      
      if (lnblocks > 0) then
        allocate(iBuff(NXB*lnblocks*3), jBuff(NYB*lnblocks*3), kBuff(NZB*lnblocks*3), dBuff(lnblocks*MDIM))
        do fc = 1, 3
          do lb = 1, lnblocks
            iBuff(NXB*lnblocks*(fc-1)+NXB*(lb-1)+1:NXB*lnblocks*(fc-1)+NXB*lb) = gr_oneBlock(lb)%firstAxisCoords(fc,gr_ilo:gr_ihi)
            jBuff(NYB*lnblocks*(fc-1)+NYB*(lb-1)+1:NYB*lnblocks*(fc-1)+NYB*lb) = gr_oneBlock(lb)%secondAxisCoords(fc,gr_jlo:gr_jhi)
            kBuff(NZB*lnblocks*(fc-1)+NZB*(lb-1)+1:NZB*lnblocks*(fc-1)+NZB*lb) = gr_oneBlock(lb)%thirdAxisCoords(fc,gr_klo:gr_khi)
            dBuff(lnblocks*(fc-1)+lb) = gr_delta(fc,lrefine(lb))
          end do
        end do
        gCoords(1:NXB,offset+1:offset+lnblocks,1:3,IAXIS) = reshape(iBuff, (/ NXB, lnblocks, 3 /))
        gCoords(1:NYB,offset+1:offset+lnblocks,1:3,JAXIS) = reshape(jBuff, (/ NYB, lnblocks, 3 /))
        gCoords(1:NZB,offset+1:offset+lnblocks,1:3,KAXIS) = reshape(kBuff, (/ NZB, lnblocks, 3 /))
        gDeltas(offset+1:offset+lnblocks,1:MDIM) = reshape(dBuff, (/ lnblocks, MDIM /))
        deallocate(iBuff, jBuff, kBuff, dBuff)
      endif

      offset = offset + localNumBlocks
    endif


    if (gr_meshMe == MASTER_PE .and. jproc /= MASTER_PE) then
    
      call MPI_Recv(localNumBlocks, 1, FLASH_INTEGER, jproc, 1, gr_meshComm, status, ierr)
      if (localNumBlocks > 0) then
        allocate(iBuff(NXB*localNumBlocks*3), jBuff(NYB*localNumBlocks*3), kBuff(NZB*localNumBlocks*3), dBuff(localNumBlocks*MDIM))
        call MPI_Recv(iBuff, NXB*localNumBlocks*3, FLASH_REAL, jproc, 1+IAXIS, gr_meshComm, status, ierr)  
        call MPI_Recv(jBuff, NYB*localNumBlocks*3, FLASH_REAL, jproc, 1+JAXIS, gr_meshComm, status, ierr)  
        call MPI_Recv(kBuff, NZB*localNumBlocks*3, FLASH_REAL, jproc, 1+KAXIS, gr_meshComm, status, ierr)  
        call MPI_Recv(dBuff,     localNumBlocks*3, FLASH_REAL, jproc, 2+MDIM , gr_meshComm, status, ierr)
        gCoords(1:NXB,offset+1:offset+localNumBlocks,1:3,IAXIS) = reshape(iBuff, (/ NXB, localNumBlocks, 3 /))
        gCoords(1:NYB,offset+1:offset+localNumBlocks,1:3,JAXIS) = reshape(jBuff, (/ NYB, localNumBlocks, 3 /))
        gCoords(1:NZB,offset+1:offset+localNumBlocks,1:3,KAXIS) = reshape(kBuff, (/ NZB, localNumBlocks, 3 /))
        gDeltas(offset+1:offset+localNumBlocks,1:MDIM) = reshape(dBuff, (/ localNumBlocks, MDIM /))
        deallocate(iBuff, jBuff, kBuff, dBuff)
      endif

      offset = offset + localNumBlocks
    endif


    if (gr_meshMe /= MASTER_PE .and. jproc == gr_meshMe) then
     
      call MPI_Send(lnblocks, 1, FLASH_INTEGER, MASTER_PE, 1, gr_meshComm, ierr)
      if (lnblocks > 0) then
        allocate(iBuff(NXB*lnblocks*3), jBuff(NYB*lnblocks*3), kBuff(NZB*lnblocks*3), dBuff(lnblocks*MDIM))
        do fc = 1, 3
          do lb = 1, lnblocks
            iBuff(NXB*lnblocks*(fc-1)+NXB*(lb-1)+1:NXB*lnblocks*(fc-1)+NXB*lb) = gr_oneBlock(lb)%firstAxisCoords(fc,gr_ilo:gr_ihi)
            jBuff(NYB*lnblocks*(fc-1)+NYB*(lb-1)+1:NYB*lnblocks*(fc-1)+NYB*lb) = gr_oneBlock(lb)%secondAxisCoords(fc,gr_jlo:gr_jhi)
            kBuff(NZB*lnblocks*(fc-1)+NZB*(lb-1)+1:NZB*lnblocks*(fc-1)+NZB*lb) = gr_oneBlock(lb)%thirdAxisCoords(fc,gr_klo:gr_khi)
            dBuff(lnblocks*(fc-1)+lb) = gr_delta(fc,lrefine(lb))
          end do
        end do
        call MPI_Send(iBuff, NXB*lnblocks*3, FLASH_REAL, MASTER_PE, 1+IAXIS, gr_meshComm, ierr)  
        call MPI_Send(jBuff, NYB*lnblocks*3, FLASH_REAL, MASTER_PE, 1+JAXIS, gr_meshComm, ierr)  
        call MPI_Send(kBuff, NZB*lnblocks*3, FLASH_REAL, MASTER_PE, 1+KAXIS, gr_meshComm, ierr)  
        call MPI_Send(dBuff,     lnblocks*3, FLASH_REAL, MASTER_PE, 2+MDIM , gr_meshComm, ierr)
        deallocate(iBuff, jBuff, kBuff, dBuff)
      endif
    
    endif


  end do


  if(gr_meshMe == MASTER_PE) then

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
        gCrdMax(a, b) = maxval(gCoords(1:d, :, b, a)) 
        gCrdMin(a, b) = minval(gCoords(1:d, :, b, a)) 
        gMtrMax(a, b) = maxval(gDeltas(:, a)) 
        gMtrMin(a, b) = minval(gDeltas(:, a)) 

        ! write dimensions
        dsetname = gCrdLbs(a, b)
        dset_dims = (/ d, gr_globalNumBlocks /)
        call h5screate_simple_f(2, dset_dims, dspc_id, error)
        call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspc_id, dset_id, error)  
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, gCoords(1:d, :, b, a), dset_dims, error)
       
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

      end do

      ! Consolodate axis "face" points
      allocate(fCoord(d+1, gr_globalNumBlocks))
      fCoord(1, :) = gCoords(1, :, LEFT_EDGE, a)
      fCoord(2:d+1, :) = gCoords(1:d, :, RIGHT_EDGE, a) 
      if (a == KAXIS .AND. d == 1) fCoord(2, :) = fCoord(1, :) + 0.000001

      ! find extreme values
      gCrdMax(a, 4) = maxval(fCoord) 
      gCrdMin(a, 4) = minval(fCoord) 
   
      ! write dimensions
      dsetname = gCrdLbs(a, 4)
      dset_dims = (/ d + 1, gr_globalNumBlocks /)
      call h5screate_simple_f(2, dset_dims, dspc_id, error)
      call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspc_id, dset_id, error)  
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, fCoord, dset_dims, error)
       
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
    deallocate(gDeltas)
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
