!!****if* source/Simulation/SimulationMain/INavierStokes/2D/bhagaWeber_mcHYPRE/Driver_evolveFlash
!!
!! NAME
!!
!!  Driver_evolveFlash
!!
!! SYNOPSIS
!!
!!  Driver_evolveFlash()
!!
!! DESCRIPTION
!!
!!  This is the main global driver for simulations that are:
!!      Spatially refined, State form, strang split
!!
!!  DOC: Driver_evolveFlash needs more explanation 
!!
!! NOTES
!!
!!  variables that begin with "dr_" like, dr_globalMe or dr_dt, dr_beginStep
!!  are stored in the data fortran module for the Driver unit, Driver_data.
!!  The "dr_" is meant to indicate that the variable belongs to the Driver Unit.
!!  all other normally named variables i, j, etc are local variables.
!!
!!
!!***

#define WRITE_TO_TECPLOT 1

subroutine Driver_evolveFlash()

  use Driver_data, ONLY: dr_globalMe, dr_nbegin,                    &
                         dr_nend, dr_dt, dr_wallClockTimeLimit, &
                         dr_tmax, dr_simTime, dr_redshift,      &
                         dr_nstep, dr_dtOld, dr_dtNew,          &
                         dr_restart, dr_elapsedWCTime
  use Driver_interface, ONLY : Driver_sourceTerms, Driver_computeDt, &
       Driver_getElapsedWCTime
  use Logfile_interface,ONLY : Logfile_stamp, Logfile_close
  use Timers_interface, ONLY : Timers_start, Timers_stop, &
                               Timers_getSummary
  use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                                Grid_getListOfBlocks, &
                                Grid_getBlkIndexLimits
  use Grid_interface, only: Grid_getCellCoords, Grid_solvePoisson, GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN

  use IO_interface,      ONLY : IO_output,IO_outputFinal

  use IO_data , ONLY : IO_checkpointFileIntervalStep, io_plotFileNumber, IO_plotFileIntervalTime, IO_plotFileIntervalStep

  implicit none

#include "constants.h"
#include "Flash.h"
 include "Flash_mpi.h"

  integer :: localNumBlocks
  integer :: blockCount
  integer :: blockList(MAXBLOCKS)
  integer :: sweepDummy
  integer :: count, firstfileflag
  real, pointer, dimension(:,:,:,:) :: solnData
  real, dimension(GRID_IHI_GC) :: xcell, xedge
  real, dimension(GRID_JHI_GC) :: ycell, yedge
  real, dimension(GRID_KHI_GC) :: zcell, zedge
  real :: pi, fact
  integer, dimension(6) :: bc_types
  real, dimension(2,6)  :: bc_values = 0.
  integer :: i,j,k
  integer ::  blockID,lb
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  logical :: endrun



  !write(*,*) "Begin Solution"
  !call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  !blockID = blockList(1)
  !call Grid_getBlkPtr(blockID, solnData, CENTER)
  !solnData(NSRC_VAR,:,:,:) = 1.0
  !write(*,*) blockID, solnData(NSRC_VAR,10,1,1)
  !call Grid_releaseBlkPtr(blockID,solnData,CENTER)



  call Logfile_stamp( 'Entering evolution loop' , '[Driver_evolveFlash]')
  call Timers_start("evolution")

  count = 0
  dr_dt = 1
  dr_nstep = 0
  firstfileflag = 0
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  call outtotecplot(dr_globalMe,dr_simtime,dr_dt,dr_nstep,count, &
                    0.0,blockList,blockCount,firstfileflag)


  call Heat_AD(blockCount, blockList, dr_simTime, dr_dt, dr_dtOld, sweepDummy)



!   pi = 3.141592654
!   fact = 1.0
!   bc_types(:) = GRID_PDE_BND_NEUMANN
!
!   do lb = 1,blockCount
!
!     blockID = blockList(lb)
!     write(*,*) "Begin Block", blockID
!
!
!
!
!     ! Get Blocks internal limits indexes:
!     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
!
!     write(*,*) "Get Data Ptr", blockID
!     call Grid_getBlkPtr(blockID,solnData,CENTER)
!
!     write(*,*) "Get Grid Data"
!     call Grid_getCellCoords(IAXIS, blockId, CENTER, .true., xcell, GRID_IHI_GC)
!     call Grid_getCellCoords(JAXIS, blockId, CENTER, .true., ycell, GRID_JHI_GC)
!     call Grid_getCellCoords(KAXIS, blockId, CENTER, .true., zcell, GRID_KHI_GC)
!
!     call Grid_getCellCoords(IAXIS, blockId, LEFT_EDGE, .true., xedge, GRID_IHI_GC)
!     call Grid_getCellCoords(JAXIS, blockId, LEFT_EDGE, .true., yedge, GRID_JHI_GC)
!     call Grid_getCellCoords(KAXIS, blockId, LEFT_EDGE, .true., zedge, GRID_KHI_GC)
!
!     write(*,*) "Zero Solution Data"
!     solnData(NSRC_VAR,:,:,:) = 0.0
!     !solnData(NFLD_VAR,:,:,:) = 0.0
!
!     write(*,*) "Create Source Functions"
!     do k=1, blkLimitsGC(HIGH,KAXIS)
!       do j=1, blkLimitsGC(HIGH,JAXIS)
!         do i=1, blkLimitsGC(HIGH,IAXIS)
!           !fact = -8.0 * pi**2 * cos(real(2.0 * pi * xcell(i))) * cos(real(2.0 * pi * ycell(j)))
!           !solnData(NSRC_VAR,i,j,k) = fact
!           write(*,*) i, j, k, solnData(NSRC_VAR,i,j,k)
!         enddo
!       enddo
!     enddo
!     fact = 1.0
!
!     write(*,*) "call the solver", solnData(NSRC_VAR, 10, 1, 1)
!     call Grid_solvePoisson(NFLD_VAR, NSRC_VAR, bc_types, bc_values, fact)
!
!     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
!
!   end do




  count = 1
  dr_dt = 1 
  dr_nstep = 1
  firstfileflag = 0
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  call outtotecplot(dr_globalMe,dr_simtime,dr_dt,dr_nstep,count, &
                    0.0,blockList,blockCount,firstfileflag)
  call IO_output(dr_simTime,dr_dt,dr_nstep+1,dr_nbegin,endRun,PLOTFILE_ONLY)

  call Timers_stop("evolution")
  call Logfile_stamp( 'Exiting evolution loop' , '[Driver_evolveFlash]')
  call Timers_getSummary(dr_nstep)
  call Logfile_stamp( "FLASH run complete.", "LOGFILE_END")
  call Logfile_close()

  return
  
end subroutine Driver_evolveFlash



