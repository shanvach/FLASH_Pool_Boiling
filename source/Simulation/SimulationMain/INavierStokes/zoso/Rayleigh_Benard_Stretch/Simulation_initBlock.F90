!!****if* source/Simulation/SimulationMain/INavierStokes/2D/LidDrivenCavity/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(in) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.
!!
!!  Reference:
!!
!! 
!! ARGUMENTS
!!
!!  blockID -          the number of the block to update
!!  myPE   -           my processor number
!!
!! 
!!
!!***

subroutine Simulation_initBlock(blockId)

  use Simulation_Data, Only : sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax, sim_init

  use Grid_interface, only  : Grid_getCellMetrics, Grid_getBlkIndexLimits, &
                              Grid_getCellCoords, &
                              Grid_getBlkCenterCoords, Grid_getBlkBoundBox, &
                              Grid_getBlkPtr, Grid_releaseBlkPtr, &
                              Grid_solvePoisson, Grid_fillGuardCells, GRID_PDE_BND_NEUMANN

  use Grid_data, only : gr_domainBC

  use Driver_Interface, only : Driver_abortFlash

  use Driver_Data, only : dr_simTime, dr_meshMe

  use ins_interface, only : ins_divergence, ins_corrector

  use HDF5

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer, intent(in) :: blockID
 
  real, pointer, dimension(:,:,:,:) :: solnData, facexData, faceyData, facezData
  
  integer :: i, j, k, nxb, nyb, nzb
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real                       :: pi
  real, dimension(GRID_IHI_GC) :: dx, dxe
  real, dimension(GRID_JHI_GC) :: dy, dye
  real, dimension(GRID_KHI_GC) :: dz, dze
  real, dimension(GRID_IHI_GC) :: xcell, xedge
  real, dimension(GRID_JHI_GC) :: ycell, yedge
  real, dimension(GRID_KHI_GC) :: zcell, zedge

  real :: xMag, xPwr, xSgn, xSym, xAmp, xFrq,       &
          yMag, yPwr, ySgn, ySym, yAmp, yFrq, rMag, &
          xFMag, xFPwr, xFSgn, yFMag, yFPwr, yFSgn, &
          xCnt, xWth, yCnt, yWth,                   &
          xSrcPol, xSrcSin, ySrcPol, ySrcSin
  integer, parameter :: seed = 86456

  ! locals necessary to read hdf5 file
  integer :: error, rank
  integer(HID_T) :: file_id, dset_id, dataspace, memspace
  integer(HSIZE_T), dimension(4) :: data_dims
  integer(HSIZE_T), dimension(4) :: count, offset, stride, block 
  real, allocatable, dimension(:,:,:,:) :: sdata
  character(len = 32) :: filename, dsetname

  ! locals necessary for solving poisson
  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)
  integer, dimension(6) :: bc_types  = GRID_PDE_BND_NEUMANN
  real, dimension(2,6)  :: bc_values = 0.0
  real                  :: poisfact  = 1.0

  integer :: ierr
  real :: mndivv,mxdivv
  real :: vecminaux(6),vecmaxaux(6),vecmin(6),vecmax(6)

  ! Point to Blocks centered variables:
  call Grid_getBlkPtr(blockID, solnData, CENTER)

  ! Point to Blocks face variables: 
  call Grid_getBlkPtr(blockID, facexData, FACEX)
  call Grid_getBlkPtr(blockID, faceyData, FACEY)

#if NDIM == 3
  call Grid_getBlkPtr(blockID, facezData, FACEZ)
#endif

  call Grid_getCellMetrics(IAXIS, blockId, CENTER, .true., dx, GRID_IHI_GC)  
  call Grid_getCellMetrics(JAXIS, blockId, CENTER, .true., dy, GRID_JHI_GC)  
  call Grid_getCellMetrics(KAXIS, blockId, CENTER, .true., dz, GRID_KHI_GC)  

  call Grid_getCellCoords(IAXIS, blockId, CENTER, .true., xcell, GRID_IHI_GC)  
  call Grid_getCellCoords(JAXIS, blockId, CENTER, .true., ycell, GRID_JHI_GC)  
  call Grid_getCellCoords(KAXIS, blockId, CENTER, .true., zcell, GRID_KHI_GC)  
  
  call Grid_getCellMetrics(IAXIS, blockId, LEFT_EDGE, .true., dxe, GRID_IHI_GC)  
  call Grid_getCellMetrics(JAXIS, blockId, LEFT_EDGE, .true., dye, GRID_JHI_GC)  
  call Grid_getCellMetrics(KAXIS, blockId, LEFT_EDGE, .true., dze, GRID_KHI_GC)  

  call Grid_getCellCoords(IAXIS, blockId, LEFT_EDGE, .true., xedge, GRID_IHI_GC)  
  call Grid_getCellCoords(JAXIS, blockId, LEFT_EDGE, .true., yedge, GRID_JHI_GC)  
  call Grid_getCellCoords(KAXIS, blockId, LEFT_EDGE, .true., zedge, GRID_KHI_GC)  
  

  ! Initialize Simulation Blocks
  select case (sim_init)


    ! Read initial conditions from file
    case (0)

      ! common parameters for hyperslab
      rank = 4
      count  = (/ 1, 1, 1, 1 /)
      offset = (/ 0, 0, 0, dr_meshMe /)
      stride = (/ 1, 1, 1, 1 /)

      ! initialize fortran hdf5 interface
      call h5open_f(error)

      ! open hdf5 file
      filename = "initBlock.h5"
      call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)
      if ( error < 0 ) call Driver_abortFlash("Please check input file for structure.")


      !!!! read in IAXIS velocity data
      ! prepare parameters
      dsetname  = "/velx"
      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, FACEX)
      nxb = blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) +1
      nyb = blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) +1
      nzb = blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) +1
      data_dims = (/ nxb, nyb, nzb, 1 /)
      block     = (/ nxb, nyb, nzb, 1 /)
      allocate(sdata(nxb, nyb, nzb, 1))

      ! read file
      call h5dopen_f(file_id, dsetname, dset_id, error)
      call h5dget_space_f(dset_id, dataspace, error)
      call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error, stride, block)
      call h5screate_simple_f(rank, data_dims, memspace, error)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, sdata, data_dims, error, memspace, dataspace)
      call h5dclose_f(dset_id, error)

      ! initialize velocity
      facexData(VELC_FACE_VAR,:,:,:) = 0.0  
      facexData(VELC_FACE_VAR, blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), & 
                               blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), & 
                               blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) & 
                               = sdata(:,:,:,1) 
      deallocate(sdata)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      
      !!!! read in JAXIS velocity data
      ! prepare parameters
      dsetname  = "/vely"
      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, FACEY)
      nxb = blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) +1
      nyb = blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) +1
      nzb = blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) +1
      data_dims = (/ nxb, nyb, nzb, 1 /)
      block     = (/ nxb, nyb, nzb, 1 /)
      allocate(sdata(nxb, nyb, nzb, 1))

      ! read file
      call h5dopen_f(file_id, dsetname, dset_id, error)
      call h5dget_space_f(dset_id, dataspace, error)
      call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error, stride, block)
      call h5screate_simple_f(rank, data_dims, memspace, error)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, sdata, data_dims, error, memspace, dataspace)
      call h5dclose_f(dset_id, error)

      ! initialize velocity
      faceyData(VELC_FACE_VAR,:,:,:) = 0.0  
      faceyData(VELC_FACE_VAR, blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), & 
                               blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), & 
                               blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) & 
                               = sdata(:,:,:,1) 
      deallocate(sdata)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if NDIM == 3
      !!!! read in KAXIS velocity data
      ! prepare parameters
      dsetname  = "/velz"
      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, FACEZ)
      nxb = blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) +1
      nyb = blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) +1
      nzb = blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) +1
      data_dims = (/ nxb, nyb, nzb, 1 /)
      block     = (/ nxb, nyb, nzb, 1 /)
      allocate(sdata(nxb, nyb, nzb, 1))

      ! read file
      call h5dopen_f(file_id, dsetname, dset_id, error)
      call h5dget_space_f(dset_id, dataspace, error)
      call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error, stride, block)
      call h5screate_simple_f(rank, data_dims, memspace, error)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, sdata, data_dims, error, memspace, dataspace)
      call h5dclose_f(dset_id, error)

      ! initialize velocity
      facezData(VELC_FACE_VAR,:,:,:) = 0.0  
      facezData(VELC_FACE_VAR, blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), & 
                               blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), & 
                               blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) & 
                               = sdata(:,:,:,1) 
      deallocate(sdata)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif

      !!!! read in temperature data
      ! prepare parameters
      dsetname  = "/temp"
      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, CENTER)
      nxb = blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) +1
      nyb = blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) +1
      nzb = blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) +1
      data_dims = (/ nxb, nyb, nzb, 1 /)
      block     = (/ nxb, nyb, nzb, 1 /)
      allocate(sdata(nxb, nyb, nzb, 1))

      ! read file
      call h5dopen_f(file_id, dsetname, dset_id, error)
      call h5dget_space_f(dset_id, dataspace, error)
      call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error, stride, block)
      call h5screate_simple_f(rank, data_dims, memspace, error)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, sdata, data_dims, error, memspace, dataspace)
      call h5dclose_f(dset_id, error)

      ! initialize velocity
      solnData(TEMP_VAR,:,:,:) = 0.0  
      solnData(TEMP_VAR, blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS), & 
                         blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS), & 
                         blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) & 
                         = sdata(:,:,:,1) 
      deallocate(sdata)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      ! close hdf5 file 
      call h5fclose_f(file_id, error)

      ! close hdf5 fortran interface
      call h5close_f(error)

      ! initialize remaining unknowns   
      solnData(PRES_VAR,:,:,:) = 0.0
      solnData(DELP_VAR,:,:,:) = 0.0
      solnData(DUST_VAR,:,:,:) = 0.0
      solnData(SDVO_VAR,:,:,:) = 0.0
      solnData(SDVC_VAR,:,:,:) = 0.0
      solnData(TVIS_VAR,:,:,:) = 0.0
      facexData(RHDS_FACE_VAR,:,:,:) = 0.0
      faceyData(RHDS_FACE_VAR,:,:,:) = 0.0
#if NDIM == 3
      facezData(RHDS_FACE_VAR,:,:,:) = 0.0
#endif

      call ins_divergence(facexData(VELC_FACE_VAR,:,:,:),&
                          faceyData(VELC_FACE_VAR,:,:,:),&
                          facezData(VELC_FACE_VAR,:,:,:),&
              blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
              blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
              blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
                          dx, dy, dz,                    &
                          solnData(SDVO_VAR,:,:,:) )


    ! Rayleigh Benard flow (Cold on Hot)
    case (1)

      solnData(PRES_VAR,:,:,:) = 0.0
      solnData(DELP_VAR,:,:,:) = 0.0
      solnData(DUST_VAR,:,:,:) = 0.0
      solnData(TVIS_VAR,:,:,:) = 0.0
      solnData(TEMP_VAR,:,:,:) = 1.0

      facexData(VELC_FACE_VAR,:,:,:) = 0.0
      faceyData(VELC_FACE_VAR,:,:,:) = 0.0
      facexData(RHDS_FACE_VAR,:,:,:) = 0.0
      faceyData(RHDS_FACE_VAR,:,:,:) = 0.0

#if NDIM == 3
      facezData(VELC_FACE_VAR,:,:,:) = 0.0
      facezData(RHDS_FACE_VAR,:,:,:) = 0.0

      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, CENTER)
      do k=1, blkLimitsGC(HIGH,KAXIS)
        do j=1, blkLimitsGC(HIGH,JAXIS)
          do i=1, blkLimitsGC(HIGH,IAXIS)
            if (zcell(k) >= 0.5) then
              solnData(TEMP_VAR,i,j,k) = 0.0
            endif
          enddo
        enddo
      enddo
#endif
#if NDIM == 2
      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, CENTER)
      do k=1, blkLimitsGC(HIGH,KAXIS)
        do j=1, blkLimitsGC(HIGH,JAXIS)
          do i=1, blkLimitsGC(HIGH,IAXIS)
            if (ycell(j) >= 0.5) then
              solnData(TEMP_VAR,i,j,k) = 0.0
            endif
          enddo
        enddo
      enddo
#endif




    ! Rayleigh Benard flow (blank)
    case (5)

      solnData(PRES_VAR,:,:,:) = 0.0
      solnData(DELP_VAR,:,:,:) = 0.0
      solnData(DUST_VAR,:,:,:) = 0.0
      solnData(TVIS_VAR,:,:,:) = 0.0
      solnData(TEMP_VAR,:,:,:) = 0.5

      facexData(VELC_FACE_VAR,:,:,:) = 0.0
      faceyData(VELC_FACE_VAR,:,:,:) = 0.0
      facexData(RHDS_FACE_VAR,:,:,:) = 0.0
      faceyData(RHDS_FACE_VAR,:,:,:) = 0.0

#if NDIM == 3
      facezData(VELC_FACE_VAR,:,:,:) = 0.0
      facezData(RHDS_FACE_VAR,:,:,:) = 0.0
#endif




    ! Rayleigh Benard flow (conv)
    case (2)

      solnData(PRES_VAR,:,:,:) = 0.0
      solnData(DELP_VAR,:,:,:) = 0.0
      solnData(DUST_VAR,:,:,:) = 0.0
      solnData(TVIS_VAR,:,:,:) = 0.0

      facexData(VELC_FACE_VAR,:,:,:) = 0.0
      faceyData(VELC_FACE_VAR,:,:,:) = 0.0
      facexData(RHDS_FACE_VAR,:,:,:) = 0.0
      faceyData(RHDS_FACE_VAR,:,:,:) = 0.0

#if NDIM == 3
      facezData(VELC_FACE_VAR,:,:,:) = 0.0
      facezData(RHDS_FACE_VAR,:,:,:) = 0.0
      
      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, CENTER)
      do k=1, blkLimitsGC(HIGH,KAXIS)
        do j=1, blkLimitsGC(HIGH,JAXIS)
          do i=1, blkLimitsGC(HIGH,IAXIS)
            solnData(TEMP_VAR,i,j,k) = -(sim_zMax/2.0)*( 2.0*(zcell(k) - 0.5) )**13.0 + sim_zMax/2.0 
          enddo
        enddo
      enddo
#endif
#if NDIM == 2
      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, CENTER)
      do k=1, blkLimitsGC(HIGH,KAXIS)
        do j=1, blkLimitsGC(HIGH,JAXIS)
          do i=1, blkLimitsGC(HIGH,IAXIS)
            solnData(TEMP_VAR,i,j,k) = -(sim_yMax/2.0)*( 2.0*(ycell(j) - 0.5) )**13.0 + sim_yMax/2.0
          enddo
        enddo
      enddo
#endif



    ! Rayleigh Benard flow (w/ Plum init)
    case (3)
                     !   RBC     IHC
                     !  ------  ------
      xMag = 1.0000  !  1.0000  1.0000  
      xPwr = 0.0000  !  0.0000  0.0000
      xSgn = 0.0000  !  0.0000  0.0000
      xSym = 0.0000  !  0.0000  0.0000
      xAmp = 0.2000  !  0.2000  0.2000
      xFrq = 2.0000  !  2.0000  2.0000

      yMag = 1.0000  !  1.0000  2.0000
      yPwr = 4.0000  !  4.0000  4.0000
      ySgn = 1.0000  !  1.0000  0.0000
      ySym = 0.0000  !  0.0000  1.0000
      yAmp = 0.0000  !  0.0000  0.0000
      yFrq = 0.0000  !  0.0000  0.0000

      rMag = 0.0025  !  0.0025  0.0025

      xFMag = 1.000  !  1.0000  1.0000
      xFPwr = 3.000  !  3.0000  3.0000
      xFSgn = 1.000  !  1.0000  1.0000
      yFMag = 1.000  !  1.0000  1.0000
      yFPwr = 2.000  !  2.0000  2.0000
      yFSgn = 0.000  !  0.0000  0.0000

      xWth =  sim_xMax - sim_xMin
      xCnt = (sim_xMax + sim_xMin) / 2.0
      yWth =  sim_yMax - sim_yMin
      yCnt = (sim_yMax + sim_yMin) / 2.0
      !call srand(seed)

      solnData(PRES_VAR,:,:,:) = 0.0
      solnData(DELP_VAR,:,:,:) = 0.0
      solnData(DUST_VAR,:,:,:) = 0.0
      solnData(TVIS_VAR,:,:,:) = 0.0

      facexData(VELC_FACE_VAR,:,:,:) = 0.0
      faceyData(VELC_FACE_VAR,:,:,:) = 0.0
      facexData(RHDS_FACE_VAR,:,:,:) = 0.0
      faceyData(RHDS_FACE_VAR,:,:,:) = 0.0

      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, CENTER)
      do k=1, blkLimitsGC(HIGH,KAXIS)
        do j=1, blkLimitsGC(HIGH,JAXIS)
          do i=1, blkLimitsGC(HIGH,IAXIS)
            
            xSrcPol = ((-1)**xSgn * (xSym * (xCnt - xcell(i) + 1e-6) / abs(xCnt - xcell(i) + 1e-6) + (1.0 - xSym)) *   &
                                    (2.0 * (xcell(i) - xCnt) / xWth)**max(0.0, 2.0 * xPwr - 1) + 1) * xMag / 2.0
            xSrcSin = xMag * xAmp * sin(2 * (2.0 * xFrq - xSgn) * 3.14 * (xcell(i) - xCnt) / xWth) *    &
                      (((-1)**xFSgn * (2.0 * (xcell(i) - xCnt) / xWth)**(2 * xFPwr) + xFSgn) * xFMag) * &
                      (((-1)**yFSgn * (2.0 * (ycell(j) - yCnt) / yWth)**(2 * yFPwr) + yFSgn) * yFMag)
            
            ySrcPol = ((-1)**ySgn * (ySym * (yCnt - ycell(j) + 1e-6) / abs(yCnt - ycell(j) + 1e-6) + (1.0 - ySym)) *   &
                                    (2.0 * (ycell(j) - yCnt) / yWth)**max(0.0, 2.0 * yPwr - 1) + 1) * yMag / 2.0
            ySrcSin = yMag * yAmp * sin(2 * (2.0 * yFrq - ySgn) * 3.14 * (ycell(j) - yCnt) / yWth) *    &
                      (((-1)**xFSgn * (2.0 * (xcell(i) - xCnt) / xWth)**(2 * xFPwr) + xFSgn) * xFMag) * &
                      (((-1)**yFSgn * (2.0 * (ycell(j) - yCnt) / yWth)**(2 * yFPwr) + yFSgn) * yFMag)
            
            !solnData(TEMP_VAR,i,j,k) = min(1.0, max(0.0, (xSrcPol + xSrcSin) * (ySrcPol + ySrcSin) + rMag * (2.0 * rand() - 1)))
            solnData(TEMP_VAR,i,j,k) = min(1.0, max(0.0, (xSrcPol + xSrcSin) * (ySrcPol + ySrcSin) + rMag * (2.0 * 0.5 - 1)))
          enddo
        enddo
      enddo
    



    ! channel flow (left to right w/ T=0)
    case (10)

      solnData(PRES_VAR,:,:,:) = 0.0
      solnData(DELP_VAR,:,:,:) = 0.0
      solnData(DUST_VAR,:,:,:) = 0.0
      solnData(TVIS_VAR,:,:,:) = 0.0
      solnData(TEMP_VAR,:,:,:) = 0.0

      facexData(VELC_FACE_VAR,:,:,:) = 1.0
      faceyData(VELC_FACE_VAR,:,:,:) = 0.0
      facexData(RHDS_FACE_VAR,:,:,:) = 0.0
      faceyData(RHDS_FACE_VAR,:,:,:) = 0.0

#if NDIM == 3
      facezData(VELC_FACE_VAR,:,:,:) = 0.0
      facezData(RHDS_FACE_VAR,:,:,:) = 0.0
#endif      



    ! channel flow (left to right w/ init T)
    case (11)

      solnData(PRES_VAR,:,:,:) = 0.0
      solnData(DELP_VAR,:,:,:) = 0.0
      solnData(DUST_VAR,:,:,:) = 0.0
      solnData(TVIS_VAR,:,:,:) = 0.0
      solnData(TEMP_VAR,:,:,:) = 0.0

      facexData(VELC_FACE_VAR,:,:,:) = 1.0
      faceyData(VELC_FACE_VAR,:,:,:) = 0.0
      facexData(RHDS_FACE_VAR,:,:,:) = 0.0
      faceyData(RHDS_FACE_VAR,:,:,:) = 0.0

#if NDIM == 3
      facezData(VELC_FACE_VAR,:,:,:) = 0.0
      facezData(RHDS_FACE_VAR,:,:,:) = 0.0
#endif      
      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, CENTER)
      do k=1, blkLimitsGC(HIGH,KAXIS)
        do j=1, blkLimitsGC(HIGH,JAXIS)
          do i=1, blkLimitsGC(HIGH,IAXIS)
            solnData(TEMP_VAR,i,j,k) = (1.0 - xcell(i))
          enddo
        enddo
      enddo




    ! channel flow (bottom to top w/ T=8)
    case (12)

      solnData(PRES_VAR,:,:,:) = 0.0
      solnData(DELP_VAR,:,:,:) = 0.0
      solnData(DUST_VAR,:,:,:) = 0.0
      solnData(TVIS_VAR,:,:,:) = 0.0
      solnData(TEMP_VAR,:,:,:) = 8.0

      facexData(VELC_FACE_VAR,:,:,:) = 0.0
      faceyData(VELC_FACE_VAR,:,:,:) = 1.0
      facexData(RHDS_FACE_VAR,:,:,:) = 0.0
      faceyData(RHDS_FACE_VAR,:,:,:) = 0.0

#if NDIM == 3
      facezData(VELC_FACE_VAR,:,:,:) = 0.0
      facezData(RHDS_FACE_VAR,:,:,:) = 0.0
#endif      



    ! channel flow (front to back w/ T=0)
    case (13)

      solnData(PRES_VAR,:,:,:) = 0.0
      solnData(DELP_VAR,:,:,:) = 0.0
      solnData(DUST_VAR,:,:,:) = 0.0
      solnData(TVIS_VAR,:,:,:) = 0.0
      solnData(TEMP_VAR,:,:,:) = 0.0

      facexData(VELC_FACE_VAR,:,:,:) = 0.0
      faceyData(VELC_FACE_VAR,:,:,:) = 0.0
      facexData(RHDS_FACE_VAR,:,:,:) = 0.0
      faceyData(RHDS_FACE_VAR,:,:,:) = 0.0

#if NDIM == 3
      facezData(VELC_FACE_VAR,:,:,:) = 1.0
      facezData(RHDS_FACE_VAR,:,:,:) = 0.0
#endif      



    ! taylor-green flow
    case (20)
      
      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, FACEX)
      do k=1, blkLimitsGC(HIGH,KAXIS)
        do j=1, blkLimitsGC(HIGH,JAXIS)
          do i=1, blkLimitsGC(HIGH,IAXIS)
            facexData(VELC_FACE_VAR,i,j,k) = -EXP(-2.0*dr_simTime)*COS(xedge(i))*SIN(ycell(j))
          enddo
        enddo
      enddo

      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, FACEY)
      do k=1, blkLimitsGC(HIGH,KAXIS)
        do j=1, blkLimitsGC(HIGH,JAXIS)
          do i=1, blkLimitsGC(HIGH,IAXIS)
            faceyData(VELC_FACE_VAR,i,j,k) = EXP(-2.0*dr_simTime)*COS(yedge(j))*SIN(xcell(i))
          enddo
        enddo
      enddo
            
      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, CENTER)
      do k=1, blkLimitsGC(HIGH,KAXIS)
        do j=1, blkLimitsGC(HIGH,JAXIS)
          do i=1, blkLimitsGC(HIGH,IAXIS)
            solnData(TEMP_VAR,i,j,k) = -EXP(-2.0*dr_simTime)*COS(xcell(i))*COS(ycell(j))/2.0 + 0.5
          enddo
        enddo
      enddo

      solnData(PRES_VAR,:,:,:) = 0.0
      solnData(DELP_VAR,:,:,:) = 0.0
      solnData(DUST_VAR,:,:,:) = 0.0
      solnData(TVIS_VAR,:,:,:) = 0.0

      facexData(RHDS_FACE_VAR,:,:,:) = 0.0
      faceyData(RHDS_FACE_VAR,:,:,:) = 0.0
    
#if NDIM == 3
      facezData(VELC_FACE_VAR,:,:,:) = 0.0
      facezData(RHDS_FACE_VAR,:,:,:) = 0.0
#endif      




    ! taylor-green flow (Modified Temp)
    case (21)
      
      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, CENTER)
      do k=1, blkLimitsGC(HIGH,KAXIS)
        do j=1, blkLimitsGC(HIGH,JAXIS)
          do i=1, blkLimitsGC(HIGH,IAXIS)
             facexData(VELC_FACE_VAR,i,j,k) = -EXP(-2.0*dr_simTime)*COS(xedge(i))*SIN(ycell(j))
          enddo
        enddo
      enddo

      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, CENTER)
      do k=1, blkLimitsGC(HIGH,KAXIS)
        do j=1, blkLimitsGC(HIGH,JAXIS)
          do i=1, blkLimitsGC(HIGH,IAXIS)
            faceyData(VELC_FACE_VAR,i,j,k) = EXP(-2.0*dr_simTime)*COS(yedge(j))*SIN(xcell(i))
          enddo
        enddo
      enddo
            
      call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC, CENTER)
      do k=1, blkLimitsGC(HIGH,KAXIS)
        do j=1, blkLimitsGC(HIGH,JAXIS)
          do i=1, blkLimitsGC(HIGH,IAXIS)
            solnData(TEMP_VAR,i,j,k) = -EXP(-2.0*dr_simTime)*SIN(xcell(i)/2.)*SIN(ycell(j)/2.)/2.0 + 0.5
          enddo
        enddo
      enddo

      solnData(PRES_VAR,:,:,:) = 0.0
      solnData(DELP_VAR,:,:,:) = 0.0
      solnData(DUST_VAR,:,:,:) = 0.0
      solnData(TVIS_VAR,:,:,:) = 0.0

      facexData(RHDS_FACE_VAR,:,:,:) = 0.0
      faceyData(RHDS_FACE_VAR,:,:,:) = 0.0

#if NDIM == 3
      facezData(VELC_FACE_VAR,:,:,:) = 0.0
      facezData(RHDS_FACE_VAR,:,:,:) = 0.0
#endif      




    
    ! null initialization  
    case default

      solnData(PRES_VAR,:,:,:) = 0.0
      solnData(DELP_VAR,:,:,:) = 0.0
      solnData(DUST_VAR,:,:,:) = 0.0
      solnData(TVIS_VAR,:,:,:) = 0.0
      solnData(TEMP_VAR,:,:,:) = 0.0

      facexData(VELC_FACE_VAR,:,:,:) = 0.0
      faceyData(VELC_FACE_VAR,:,:,:) = 0.0
      facexData(RHDS_FACE_VAR,:,:,:) = 0.0
      faceyData(RHDS_FACE_VAR,:,:,:) = 0.0

#if NDIM == 3
      facezData(VELC_FACE_VAR,:,:,:) = 0.0
      facezData(RHDS_FACE_VAR,:,:,:) = 0.0
#endif      




  
  end select

  call Grid_releaseBlkPtr(blockID,solnData,CENTER)
  call Grid_releaseBlkPtr(blockID,facexData,FACEX)
  call Grid_releaseBlkPtr(blockID,faceyData,FACEY)

#if NDIM == 3
  call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#endif


  return

111    format (i4,3x,i4)
112    format (3(3x,e12.4))

end subroutine Simulation_initBlock
