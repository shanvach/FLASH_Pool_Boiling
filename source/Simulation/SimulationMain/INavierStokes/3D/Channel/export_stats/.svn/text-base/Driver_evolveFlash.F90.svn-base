!!****if* source/Simulation/SimulationMain/INavierStokes/2D/LidDrivenCavity/Driver_evolveFlash
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
!!$#define DEDUG_ALL

!!#define WRITE_TO_TECPLOT 1

#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif

subroutine Driver_evolveFlash()

  use Driver_data, ONLY: dr_globalMe, dr_nbegin,       &
                         dr_nend, dr_dt, dr_wallClockTimeLimit, &
                         dr_tmax, dr_simTime, dr_redshift,      &
                         dr_nstep, dr_dtOld, dr_dtNew,          &
                         dr_restart, dr_elapsedWCTime
  use IncompNS_interface, ONLY : IncompNS
  use Driver_interface, ONLY : Driver_sourceTerms, Driver_computeDt, &
       Driver_getElapsedWCTime
  use Logfile_interface,ONLY : Logfile_stamp, Logfile_close
  use Timers_interface, ONLY : Timers_start, Timers_stop, &
                               Timers_getSummary
  use Particles_interface, ONLY : Particles_advance, Particles_dump

  use Grid_interface, ONLY : Grid_getListOfBlocks,   &
                             Grid_getDeltas,         &
                             Grid_getBlkBC,          &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkIndexLimits, &
                             Grid_fillGuardCells,    &
                             Grid_putFluxData,       &
                             Grid_getFluxData,       &
                             Grid_conserveFluxes,    &
                             Grid_conserveField,     &
                             Grid_solvePoisson, Grid_getBlkBoundBox, Grid_getBlkCenterCoords

  use gr_interface, ONLY : gr_findMean


  use Gravity_interface, ONLY :  Gravity_potentialListOfBlocks
  use IO_interface,      ONLY : IO_output,IO_outputFinal

  use Profiler_interface, ONLY : Profiler_start, Profiler_stop

#ifdef WRITE_TO_TECPLOT
  use IO_data , ONLY : io_plotFileNumber, IO_plotFileIntervalStep, IO_plotFileIntervalTime
#endif

  use ins_interface, only  :  ins_velomg2center, &
                              ins_divergence, &
                              ins_fluxfix,    &
                              ins_fluxfix_p,  &
                              ins_corrector

  use IncompNS_data, only : ins_cflflg


  use Grid_data, only : gr_meshMe,gr_meshComm

  use IncompNSstats_interface, only : IncompNSstats
 
  use instats_interface, only : instats_ioexport

  implicit none

#include "constants.h"
#include "Flash.h"
 include "Flash_mpi.h"

  integer   :: localNumBlocks

  integer :: blockCount
  integer :: blockList(MAXBLOCKS)
  integer :: sweepDummy

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

  integer :: lb,blockID,i,j,jj,k

  
  ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(3,2) :: strBuff
  character(len=15) :: numToStr
  character(20) :: filename

  logical :: endRun

  logical :: tecplot_flg
  integer count, firstfileflag
  
  endRun = .false.

  call Logfile_stamp( 'Entering evolution loop' , '[Driver_evolveFlash]')
  call Profiler_start("FLASH_evolution")
  call Timers_start("evolution")

  ! Initial Timestep:
  ! backup needed old
  dr_dtOld = dr_dt

  ! calculate new
  call Driver_computeDt(dr_nbegin,  dr_nstep,      &
                        dr_simTime, dr_dtOld, dr_dtNew)
  ! store new
  dr_dt = dr_dtNew

  if (dr_globalMe == MASTER_PE) write(*,*) 'dr_dt ===',dr_dt

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

#ifdef WRITE_TO_TECPLOT
  count = io_plotFileNumber-1
  if (.not. dr_restart) then
  firstfileflag = 0
  dr_nstep = 0
  count    = 0
  call outtotecplot(dr_globalMe,dr_simtime,dr_dt,dr_nstep,count, &
                    0.0,blockList,blockCount,firstfileflag)
  endif
  firstfileflag = 1
#endif

  dr_nstep = dr_nbegin

  ! Call stats export routine:
  call instats_ioexport(.true.) 

  ! Write binaries for each block:
  write(filename, '("IOData/Restrt.", i5.5)') gr_meshMe
  open(33,file=trim(filename),form='unformatted')

  ! Get List of leaf Blocks:
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  ! Loop on Blocks
  do lb = 1,blockCount
     blockID = blockList(lb)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,facezData,FACEZ)


     ! U velocity:
     do k=GRID_KLO,GRID_KHI
       do j=GRID_JLO,GRID_JHI
         do i=GRID_JLO,GRID_IHI+1

           write(33) facexData(VELC_FACE_VAR,i,j,k)

         enddo
       enddo
     enddo

     ! V velocity:
     do k=GRID_KLO,GRID_KHI
       do j=GRID_JLO,GRID_JHI+1
         do i=GRID_JLO,GRID_IHI

           write(33) faceyData(VELC_FACE_VAR,i,j,k)

         enddo
       enddo
     enddo

     ! W velocity:
     do k=GRID_KLO,GRID_KHI+1
       do j=GRID_JLO,GRID_JHI
         do i=GRID_JLO,GRID_IHI

           write(33) facezData(VELC_FACE_VAR,i,j,k)

         enddo
       enddo
     enddo

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  enddo

  close(33)

  call Timers_stop("evolution")
  call Profiler_stop("FLASH_evolution")
  call Logfile_stamp( 'Exiting evolution loop' , '[Driver_evolveFlash]')
  call Timers_getSummary(dr_nstep)
  call Logfile_stamp( "FLASH run complete.", "LOGFILE_END")
  call Logfile_close()

  return
  
end subroutine Driver_evolveFlash



