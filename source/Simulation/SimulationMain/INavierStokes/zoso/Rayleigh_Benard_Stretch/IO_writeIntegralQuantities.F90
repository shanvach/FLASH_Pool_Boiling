!!****if* source/IO/IOMain/IO_writeIntegralQuantities
!!
!!
!!  NAME
!!    IO_writeIntegralQuantities
!!
!!  SYNOPSIS
!!    call IO_writeIntegralQuantities(integer(in) :: isFirst,
!!                                    real(in)    :: simTime)
!!
!!  DESCRIPTION
!!
!!   Compute the values of integral quantities (eg. total energy)
!!   and write them to an ASCII file.  If this is the initial step,
!!   create the file and write a header to it before writing the data.
!!
!!  ARGUMENTS
!!    
!!   isFirst - if 1 then write header info plus data, otherwise just write data
!!   simTime - simulation time
!!
!!
!!***

subroutine IO_writeIntegralQuantities ( isFirst, simTime)

  use IO_data, ONLY : io_globalMe, io_restart, io_statsFileName, io_globalComm
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getDomainBoundBox,                &
                             Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_getCellMetrics, &
                             Grid_releaseBlkPtr
  use Grid_data, ONLY : gr_axisNumProcs, gr_axisMe
  use Heat_AD_data, ONLY : ht_invsqrtRaPr
  use IncompNS_data, ONLY : ins_invsqrtRa_Pr

  implicit none

#include "Flash_mpi.h"
#include "constants.h"
#include "Flash.h"
    
  real, intent(in) :: simTime
  integer, intent(in) :: isFirst
  
  integer :: funit = 99
  integer :: error
  character (len=MAX_STRING_LENGTH), save :: fname 
  integer :: lb, blockCount, blockID
  integer :: blockList(MAXBLOCKS)
  integer :: blkLimits(HIGH, MDIM), blkLimitsGC(HIGH, MDIM)
  integer, parameter ::  nGlobalSum = 8  ! Number of globally-summed quantities
  integer, parameter ::  nScreenSum = 10 ! Number of screen output quantities
  real :: gsum(nGlobalSum) ! Global summed quantities
  real :: lsum(nGlobalSum) ! Local  summed quantities
  real :: psum(nScreenSum) ! Screen output quantities
  integer :: i, j, k, nI, nJ, nK
  real :: volume, area, thermal, kinetic, masserr, convect, conduct, Nu 
  real :: qpp_lwr, qpp_upr, heatFlux, viscous, bodyWork
  real :: thermalError, thermalErrorRel, kineticError, kineticErrorRel, heatTrans
  real :: sqrtRaPr, invsqrtRaPr, invsqrtRa_Pr, bndBox(2,MDIM) 
  real, dimension(:,:,:,:), pointer :: solnData, facexData, faceyData, facezData
  real, dimension(GRID_IHI_GC,3) :: iMetrics
  real, dimension(GRID_JHI_GC,3) :: jMetrics
  real, dimension(GRID_KHI_GC,3) :: kMetrics
  integer :: ioStat
  real, save :: intHeatTrans, oldHeatTrans, intBodyWork, oldBodyWork, intViscous, oldViscous
  real, save :: oldSimTime, thermalInitial, kineticInitial, thermalScale, kineticScale

  if (io_globalMe == MASTER_PE .and. isFirst == 1) then
    intHeatTrans = 0.0
    intBodyWork  = 0.0
    intViscous   = 0.0
    oldSimTime   = 0.0
    oldHeatTrans = 0.0
    oldBodyWork  = 0.0
    oldViscous   = 0.0
    thermalScale = 1.0
    kineticScale = 1.0
  endif

  call Grid_getListOfBlocks(LEAF, blockList, blockCount)
  
  do lb = 1, blockCount
    blockID = blockList(lb)
     
    !get the index limits of the block
    call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
     
    !get the block metrics to calulate relavent volumes
    call Grid_getCellMetrics(IAXIS,blockID,LEFT_EDGE, .true.,iMetrics(:,LEFT_EDGE), GRID_IHI_GC)
    call Grid_getCellMetrics(IAXIS,blockID,CENTER,    .true.,iMetrics(:,CENTER),    GRID_IHI_GC)
    call Grid_getCellMetrics(IAXIS,blockID,RIGHT_EDGE,.true.,iMetrics(:,RIGHT_EDGE),GRID_IHI_GC)
    call Grid_getCellMetrics(JAXIS,blockID,LEFT_EDGE, .true.,jMetrics(:,LEFT_EDGE), GRID_JHI_GC)
    call Grid_getCellMetrics(JAXIS,blockID,CENTER,    .true.,jMetrics(:,CENTER),    GRID_JHI_GC)
    call Grid_getCellMetrics(JAXIS,blockID,RIGHT_EDGE,.true.,jMetrics(:,RIGHT_EDGE),GRID_JHI_GC)
    call Grid_getCellMetrics(KAXIS,blockID,LEFT_EDGE, .true.,kMetrics(:,LEFT_EDGE), GRID_KHI_GC)
    call Grid_getCellMetrics(KAXIS,blockID,CENTER,    .true.,kMetrics(:,CENTER),    GRID_KHI_GC)
    call Grid_getCellMetrics(KAXIS,blockID,RIGHT_EDGE,.true.,kMetrics(:,RIGHT_EDGE),GRID_KHI_GC)

    ! get a pointer to the current block of data
    call Grid_getBlkPtr(blockID, solnData)
    call Grid_getBlkPtr(blockID,facexData,FACEX)
    call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
    call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif

    ! Initially zero required quantities
    psum(:) = 0.0
    gsum(:) = 0.0
    lsum(:) = 0.0
    convect = 0.0
    conduct = 0.0
    qpp_lwr = 0.0
    qpp_upr = 0.0

    ! Save non-dimentional quantities
    sqrtRaPr = 1.0 / ht_invsqrtRaPr
    invsqrtRaPr = ht_invsqrtRaPr
    invsqrtRa_Pr = ins_invsqrtRa_Pr

    nI = gr_axisNumProcs(IAXIS) * (blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1)
    nJ = gr_axisNumProcs(JAXIS) * (blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1)
    nK = gr_axisNumProcs(KAXIS) * (blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1)

    ! Sum contributions from the indicated blkLimits of cells.
    do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
      do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
          ! Differential Volume
          volume = 1.0 / (iMetrics(i,CENTER) * jMetrics(j,CENTER) * kMetrics(k,CENTER))

          ! Differential Area
#if NDIM == 3
          area = 1.0 / (iMetrics(i,CENTER) * jMetrics(j,CENTER))
#else
          area = 1.0 / iMetrics(i,CENTER) 
#endif

          ! Thermal Energy 
          thermal = solnData(TEMP_VAR,i,j,k) * volume

          ! Kinetic Energy 
          kinetic = ( 0.5 * (facexData(VELC_FACE_VAR,i+1,j,k)**2 + facexData(VELC_FACE_VAR,i,j,k)**2) )
          kinetic = ( 0.5 * (faceyData(VELC_FACE_VAR,i,j+1,k)**2 + faceyData(VELC_FACE_VAR,i,j,k)**2) ) + kinetic
#if NDIM == 3
          kinetic = ( 0.5 * (facezData(VELC_FACE_VAR,i,j,k+1) + facezData(VELC_FACE_VAR,i,j,k)) )**2 + kinetic
#endif
          kinetic = 0.5 * kinetic * volume

          ! Volume weighted L1 norm of Divergence (Mass Error)
          masserr = (facexData(VELC_FACE_VAR,i+1,j,k) - facexData(VELC_FACE_VAR,i,j,k)) * iMetrics(i,CENTER)
          masserr = (faceyData(VELC_FACE_VAR,i,j+1,k) - faceyData(VELC_FACE_VAR,i,j,k)) * jMetrics(j,CENTER) + masserr
#if NDIM == 3
          masserr = (facezData(VELC_FACE_VAR,i,j,k+1) - facezData(VELC_FACE_VAR,i,j,k)) * kMetrics(k,CENTER) + masserr
#endif
          masserr = abs(masserr) * volume

          ! Planer Nusselt Number Summation
          convect = sqrtRaPr * 0.5 * (faceyData(VELC_FACE_VAR,i,j+1,k) + faceyData(VELC_FACE_VAR,i,j,k)) * solnData(TEMP_VAR,i,j,k) / nI + convect
          conduct = 0.5 * (solnData(TEMP_VAR,i,j+1,k) - solnData(TEMP_VAR,i,j-1,k)) * jMetrics(j,CENTER) / nI + conduct

          ! Calculate Surface Heat Flux
          qpp_lwr = 0.0
          qpp_upr = 0.0
#if NDIM == 3
          if(gr_axisMe(KAXIS) == 0 .and. k == blkLimits(LOW,KAXIS)) then
            qpp_lwr = -invsqrtRaPr * (solnData(TEMP_VAR,i,j,k) - solnData(TEMP_VAR,i,j,k-1)) * kMetrics(k,LEFT_EDGE) * area
          endif
          if(gr_axisMe(KAXIS) == gr_axisNumProcs(KAXIS)-1 .and. k == blkLimits(HIGH,KAXIS)) then
            qpp_upr = invsqrtRaPr * (solnData(TEMP_VAR,i,j,k+1) - solnData(TEMP_VAR,i,j,k)) * kMetrics(k,RIGHT_EDGE) * area
          endif
#else
          if(gr_axisMe(JAXIS) == 0 .and. j == blkLimits(LOW,JAXIS)) then
            qpp_lwr = -invsqrtRaPr * (solnData(TEMP_VAR,i,j,k) - solnData(TEMP_VAR,i,j-1,k)) * jMetrics(j,LEFT_EDGE) * area 
          endif
          if(gr_axisMe(JAXIS) == gr_axisNumProcs(JAXIS)-1 .and. j == blkLimits(HIGH,JAXIS)) then
            qpp_upr = invsqrtRaPr * (solnData(TEMP_VAR,i,j+1,k) - solnData(TEMP_VAR,i,j,k)) * jMetrics(j,RIGHT_EDGE) * area
          endif
#endif

          ! Calculate Body Force Work (u dot g_hat T)
#if NDIM == 3
          bodyWork = 0.5 * (facezData(VELC_FACE_VAR,i,j,k+1) + facezData(VELC_FACE_VAR,i,j,k)) * solnData(TEMP_VAR,i,j,k)
#else
          bodyWork = 0.5 * (faceyData(VELC_FACE_VAR,i,j+1,k) + faceyData(VELC_FACE_VAR,i,j,k)) * solnData(TEMP_VAR,i,j,k)
#endif
          bodyWork = bodyWork * volume

          ! Calculate Viscous Dissipation Work in two or three dimensions
          viscous = 0.25 * (facexData(VELC_FACE_VAR,i+1,j,k) + facexData(VELC_FACE_VAR,i,j,k)) * ( (                                                      &
                           (facexData(VELC_FACE_VAR,i+2,j,k) + facexData(VELC_FACE_VAR,i+1,j,k)) * iMetrics(i,RIGHT_EDGE) -                               &
                           (facexData(VELC_FACE_VAR,i+1,j,k) + facexData(VELC_FACE_VAR,i,j,k)) * (iMetrics(i,LEFT_EDGE) + iMetrics(i,RIGHT_EDGE)) + & 
                           (facexData(VELC_FACE_VAR,i,j,k) + facexData(VELC_FACE_VAR,i-1,j,k)) * iMetrics(i,LEFT_EDGE) ) * iMetrics(i,CENTER) + (         &
                           (facexData(VELC_FACE_VAR,i+1,j+1,k) + facexData(VELC_FACE_VAR,i,j+1,k)) * jMetrics(j,RIGHT_EDGE) -                             &
                           (facexData(VELC_FACE_VAR,i+1,j,k) + facexData(VELC_FACE_VAR,i,j,k)) * (jMetrics(j,LEFT_EDGE) + jMetrics(j,RIGHT_EDGE)) + & 
                           (facexData(VELC_FACE_VAR,i+1,j-1,k) + facexData(VELC_FACE_VAR,i,j-1,k)) * jMetrics(j,LEFT_EDGE) ) * jMetrics(j,CENTER) ) +     &
                    0.25 * (faceyData(VELC_FACE_VAR,i,j+1,k) + faceyData(VELC_FACE_VAR,i,j,k)) * ( (                                                      &
                           (faceyData(VELC_FACE_VAR,i+1,j+1,k) + faceyData(VELC_FACE_VAR,i+1,j,k)) * iMetrics(i,RIGHT_EDGE) -                             &
                           (faceyData(VELC_FACE_VAR,i,j+1,k) + faceyData(VELC_FACE_VAR,i,j,k)) * (iMetrics(i,LEFT_EDGE) + iMetrics(i,RIGHT_EDGE)) + & 
                           (faceyData(VELC_FACE_VAR,i-1,j+1,k) + faceyData(VELC_FACE_VAR,i-1,j,k)) * iMetrics(i,LEFT_EDGE) ) * iMetrics(i,CENTER) + (     &
                           (faceyData(VELC_FACE_VAR,i,j+2,k) + faceyData(VELC_FACE_VAR,i,j+1,k)) * jMetrics(j,RIGHT_EDGE) -                               &
                           (faceyData(VELC_FACE_VAR,i,j+1,k) + faceyData(VELC_FACE_VAR,i,j,k)) * (jMetrics(j,LEFT_EDGE) + jMetrics(j,RIGHT_EDGE)) + & 
                           (faceyData(VELC_FACE_VAR,i,j,k) + faceyData(VELC_FACE_VAR,i,j-1,k)) * jMetrics(j,LEFT_EDGE) ) * jMetrics(j,CENTER) )
#if NDIM == 3
          viscous = 0.25 * (facexData(VELC_FACE_VAR,i+1,j,k) + facexData(VELC_FACE_VAR,i,j,k)) * ( (                                                      &
                           (facexData(VELC_FACE_VAR,i+1,j,k+1) + facexData(VELC_FACE_VAR,i,j,k+1)) * kMetrics(k,RIGHT_EDGE) -                             &
                           (facexData(VELC_FACE_VAR,i+1,j,k) + facexData(VELC_FACE_VAR,i,j,k)) * (kMetrics(k,LEFT_EDGE) + kMetrics(k,RIGHT_EDGE)) + & 
                           (facexData(VELC_FACE_VAR,i+1,j,k-1) + facexData(VELC_FACE_VAR,i,j,k-1)) * kMetrics(k,LEFT_EDGE) ) * kMetrics(k,CENTER) ) +     &
                    0.25 * (faceyData(VELC_FACE_VAR,i,j+1,k) + faceyData(VELC_FACE_VAR,i,j,k)) * ( (                                                      &
                           (faceyData(VELC_FACE_VAR,i,j+1,k+1) + faceyData(VELC_FACE_VAR,i,j,k+1)) * kMetrics(k,RIGHT_EDGE) -                             &
                           (faceyData(VELC_FACE_VAR,i,j+1,k) + faceyData(VELC_FACE_VAR,i,j,k)) * (kMetrics(k,LEFT_EDGE) + kMetrics(k,RIGHT_EDGE)) + & 
                           (faceyData(VELC_FACE_VAR,i,j+1,k-1) + faceyData(VELC_FACE_VAR,i,j,k-1)) * kMetrics(k,LEFT_EDGE) ) * kMetrics(k,CENTER) ) +     &
                    0.25 * (facezData(VELC_FACE_VAR,i,j,k+1) + facezData(VELC_FACE_VAR,i,j,k)) * ( (                                                      &
                           (facezData(VELC_FACE_VAR,i+1,j,k+1) + facezData(VELC_FACE_VAR,i+1,j,k)) * iMetrics(i,RIGHT_EDGE) -                             &
                           (facezData(VELC_FACE_VAR,i,j,k+1) + facezData(VELC_FACE_VAR,i,j,k)) * (iMetrics(i,LEFT_EDGE) + iMetrics(i,RIGHT_EDGE)) + & 
                           (facezData(VELC_FACE_VAR,i-1,j,k+1) + facezData(VELC_FACE_VAR,i-1,j,k)) * iMetrics(i,LEFT_EDGE) ) * iMetrics(i,CENTER) + (     &
                           (facezData(VELC_FACE_VAR,i,j+1,k+1) + facezData(VELC_FACE_VAR,i,j+1,k)) * jMetrics(j,RIGHT_EDGE) -                             &
                           (facezData(VELC_FACE_VAR,i,j,k+1) + facezData(VELC_FACE_VAR,i,j,k)) * (jMetrics(j,LEFT_EDGE) + jMetrics(j,RIGHT_EDGE)) + & 
                           (facezData(VELC_FACE_VAR,i,j-1,k+1) + facezData(VELC_FACE_VAR,i,j-1,k)) * jMetrics(j,LEFT_EDGE) ) * jMetrics(j,CENTER) + (     &
                           (facezData(VELC_FACE_VAR,i,j,k+2) + facezData(VELC_FACE_VAR,i,j,k+1)) * kMetrics(k,RIGHT_EDGE) -                               &
                           (facezData(VELC_FACE_VAR,i,j,k+1) + facezData(VELC_FACE_VAR,i,j,k)) * (kMetrics(k,LEFT_EDGE) + kMetrics(k,RIGHT_EDGE)) + & 
                           (facezData(VELC_FACE_VAR,i,j,k) + facezData(VELC_FACE_VAR,i,j,k-1)) * kMetrics(k,LEFT_EDGE) ) * kMetrics(k,CENTER) ) + viscous
#endif
          viscous = invsqrtRa_Pr * viscous * volume

          ! Save quantities integrated over volume
          lsum(1) = lsum(1) + thermal
          lsum(2) = lsum(2) + kinetic
          lsum(3) = lsum(3) + masserr
          lsum(4) = lsum(4) + qpp_lwr
          lsum(5) = lsum(5) + qpp_upr
          lsum(6) = lsum(6) + bodyWork
          lsum(7) = lsum(7) + viscous

        end do
       
#if NDIM == 2
        ! 2D Average Planer Nusselt Number
        Nu = (convect - conduct) / nJ
        convect = 0.0
        conduct = 0.0
        lsum(8) = lsum(8) + Nu
#endif

      end do

#if NDIM == 3
        ! 3D Average Planer Nusselt Number
        Nu = (convect - conduct) / nJ
        convect = 0.0
        conduct = 0.0
        lsum(8) = lsum(8) + Nu
#endif

    end do

    call Grid_releaseBlkPtr(blockID,solnData,CENTER)
    call Grid_releaseBlkPtr(blockID,facexData,FACEX)
    call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
    call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

  enddo
  
  call MPI_Reduce(lsum, gsum, nGlobalSum, FLASH_REAL, MPI_SUM, MASTER_PE, io_globalComm, error)
  
  if (io_globalMe  == MASTER_PE) then
    
    ! Calculate Integral Geometric Quantities
    call Grid_getDomainBoundBox(bndBox)
    volume = (bndBox(HIGH,IAXIS) - bndBox(LOW,IAXIS)) * (bndBox(HIGH,JAXIS) - bndBox(LOW,JAXIS))
    area = bndBox(HIGH,IAXIS) - bndBox(LOW,IAXIS)
#if NDIM == 3
    volume = volume * (bndBox(HIGH,KAXIS) - bndBox(LOW,KAXIS))
    area = area * (bndBox(HIGH,JAXIS) - bndBox(LOW,JAXIS))
#endif

    ! Normalize desired quantities
    gsum(1) = gsum(1) / volume  ! Thermal Energy [energy / volume]
    gsum(2) = gsum(2) / volume  ! Kinetic Energy [energy / volume]
    gsum(3) = gsum(3) / volume  ! Mass Error
    gsum(4) = gsum(4) / area    ! Heat Flux [energy / area time]
    gsum(5) = gsum(5) / area    ! Heat Flux [energy / area time]
    gsum(6) = gsum(6) / volume  ! Body Work [energy / volume]
    gsum(7) = gsum(7) / volume  ! viscous Work [energy / volume]

    ! Calculate Integral Heat Transfer 
    heatTrans = (gsum(4) + gsum(5)) * area / volume
    intHeatTrans = intHeatTrans + (simTime - oldSimTime) * (heatTrans + oldHeatTrans) / 2.0
    oldHeatTrans = heatTrans

    ! Calculate Integral Body Work
    intBodyWork = intBodyWork + (simTime - oldSimTime) * (gsum(6) + oldBodyWork) / 2.0
    oldBodyWork = gsum(6)

    ! Calculate Integral Viscous Work
    intViscous = intViscous + (simTime - oldSimTime) * (gsum(7) + oldViscous) / 2.0
    oldViscous = gsum(7)

    ! Save Last Simulation Time for Next Integrals
    oldSimTime = simTime

    ! Calculate Relative Error in Thermal Energy
    if (isFirst == 1) then
            thermalInitial = gsum(1)
            if (abs(thermalInitial) > 1E-4) thermalScale = thermalInitial
    endif
    thermalError = gsum(1) - intHeatTrans - thermalInitial
    thermalErrorRel = 100.0 * abs(thermalError) / thermalScale 

    ! Calculate Relative Error in Kinetic Energy
    if (isFirst == 1) then
            kineticInitial = gsum(2)
            if (abs(kineticInitial) > 1E-4) kineticScale = kineticScale
    endif
    kineticError = gsum(2) - intBodyWork -intViscous - kineticInitial
    kineticErrorRel = 100.0 * abs(kineticError) / kineticScale 

    ! Gather Quantities for Screen Output
    psum(1)  = simTime         ! Simulation Time
    psum(2)  = gsum(1)         ! Thermal Energy
    psum(3)  = gsum(2)         ! Kinetic Energy
    psum(4)  = gsum(8)         ! Nusselt Number
    psum(5)  = thermalErrorRel ! Thermal Energy Error (Relative)
    psum(6)  = kineticErrorRel ! kinetic Energy Error (Relative)
    psum(7)  = gsum(3)         ! Mass Error (Volume wt L1 Div U)
    psum(8)  = intHeatTrans    ! Integral Heat Transfer
    psum(9)  = oldBodyWork     ! Integral Body Work
    psum(10) = oldViscous      ! Integral viscous Work

    ! create the file from scratch if it is a not a restart simulation, 
    ! otherwise append to the end of the file
    ioStat = 0
    open(funit, file=trim(io_statsFileName), position='APPEND', status='OLD', iostat=ioStat)
    if(ioStat .NE. 0) then
      open(funit, file=trim(io_statsFileName), position='APPEND')
    endif
    
    ! Write header to file
    if (isFirst .EQ. 1 .AND. (.NOT. io_restart .or. ioStat .NE. 0)) then
      write (funit, 10)               &
        'Simulation time           ', &
        'Energy -- Thermal         ', &
        'Energy -- Kinetic         ', &
        'Nusselt (Nu) Number       ', &
        'Thermal Error (Rel) [%]   ', &
        'Kinetic Error (Rel) [%]   ', &
        'Mass Error (Vwt L1 divU)  ', &
        'Int Heat Trans (Int Q)    ', &
        'Int Body Work (Int uT)    ', &
        'Int Viscous Work          '
10      format (2x,50(a25, :, 1X))
    else if(isFirst .EQ. 1) then
      write (funit, 11) 
11      format('# simulation restarted')
    endif
     
    ! Write the global sums to the file.
    write (funit, 12) psum 
12  format (1x, 50(es25.18, :, 1x))
 
    close (funit)          ! Close the file.
     
  endif
  
  call MPI_Barrier(io_globalComm, error)
  
  return
end subroutine IO_writeIntegralQuantities



