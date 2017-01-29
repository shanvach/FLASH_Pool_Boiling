!!****if* source/IO/IOMain/MPIHybrid/IO_init
!!
!! NAME
!!  IO_init
!!
!! SYNOPSIS
!!
!!  call IO_init()
!!
!! DESCRIPTION
!!
!!  Initialize the data module for the IO unit using hybrid parallel
!!  sequential MPI-IO implementation
!!
!! ARGUMENTS
!!
!!
!!***



subroutine IO_init()
  
  use IO_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
    RuntimeParameters_mapStrToInt, RuntimeParameters_getNumReal, &
    RuntimeParameters_getNumInt, RuntimeParameters_getNumStr, RuntimeParameters_getNumLog
  use Grid_interface, ONLY : Grid_getBlkIndexLimits
  use IO_interface, ONLY : IO_readCheckpoint
  use IO_data, ONLY : io_globalMe, io_glabalNumProcs
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"


  integer :: i, count, ierr, siz, fh
  
  integer :: gsize,lsize,start
  integer :: color, key,arrayDim
  integer :: blockID=1
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: i

  io_lastWallClockCheckpoint = MPI_Wtime()

  io_CPUSeconds = 0.0


  call Driver_getMype(GLOBAL_COMM, io_globalMe)

  call Driver_getNumProcs(GLOBAL_COMM, io_globalNumProcs)

  call RuntimeParameters_get("procGroup",io_procGroup)
  call RuntimeParameters_get('plotFileNumber', io_plotFileNumber)

  call RuntimeParameters_get('checkpointFileNumber', io_checkpointFileNumber)

  call RuntimeParameters_get('tinitial', io_tinitial)
  call RuntimeParameters_get('plotfileIntervalTime',    io_plotfileIntervalTime)
  call RuntimeParameters_get('checkpointFileIntervalTime',   io_checkpointFileIntervalTime)
  call RuntimeParameters_get('checkpointFileIntervalStep',   io_checkpointFileIntervalStep)

  call RuntimeParameters_get('zinitial', io_zinitial)
  call RuntimeParameters_get('zplot',    io_zplot)
  call RuntimeParameters_get('zrstrt',   io_zrstrt)
  
  call RuntimeParameters_get('restart',  io_restart)
  call RuntimeParameters_get('rolling_checkpoint', io_rollingCheckpoint)
  call RuntimeParameters_get('wall_clock_checkpoint', io_wallClockCheckpoint)

  !DEV: what to do with memoryStatFreq not yet implemented
  call RuntimeParameters_get('memory_stat_freq', io_memoryStatFreq)

  call RuntimeParameters_get('wr_integrals_freq', io_integralFreq)
  call RuntimeParameters_get('stats_file', io_statsFileName)

  ! get the runtime parameters
  
  call RuntimeParameters_get('basenm', io_baseName)
  call RuntimeParameters_get("output_directory", io_outputDir)
  call RuntimeParameters_get('geometry', io_geometry)


#ifdef FLASH_GRID_UG
  call RuntimeParameters_get('iprocs', io_iprocs)
  call RuntimeParameters_get('jprocs', io_jprocs)
  call RuntimeParameters_get('kprocs', io_kprocs)
#endif

  ! get all the IO plot vars.  These are the varibles that will be written 
  ! to plotfiles
  call RuntimeParameters_get("plot_var_1", io_plotVarStr(1))
  call RuntimeParameters_get('plot_var_2', io_plotVarStr(2))
  call RuntimeParameters_get('plot_var_3', io_plotVarStr(3))
  call RuntimeParameters_get('plot_var_4', io_plotVarStr(4))
  call RuntimeParameters_get('plot_var_5', io_plotVarStr(5))
  call RuntimeParameters_get('plot_var_6', io_plotVarStr(6))
  call RuntimeParameters_get('plot_var_7', io_plotVarStr(7))
  call RuntimeParameters_get('plot_var_8', io_plotVarStr(8))
  call RuntimeParameters_get('plot_var_9', io_plotVarStr(9))
  call RuntimeParameters_get('plot_var_10', io_plotVarStr(10))
  call RuntimeParameters_get('plot_var_11', io_plotVarStr(11))
  call RuntimeParameters_get('plot_var_12', io_plotVarStr(11))


  call RuntimeParameters_get('bytePack', io_bytePack)

#ifdef FLASH_PARTICLES

  call RuntimeParameters_get('particleFileIntervalTime', io_particleFileIntevalTime)
  call RuntimeParameters_get('particleFileNumber', io_particleFileNumber)
#endif

  !translate the plotvars from strings to integers
  !and calculate the number of them
  io_nPlotVars = 0
  do i=UNK_VARS_BEGIN,UNK_VARS_END
     call RuntimeParameters_mapStrToInt(trim(io_plotVarStr(i)), io_plotVar(i))
     if(io_plotVar(i) /= NONEXISTENT) then
        io_nPlotVars = io_nPlotVars + 1
     end if
  end do

  

  if(.not. io_restart) then

     call RuntimeParameters_getNumReal(io_numRealParms)
     call RuntimeParameters_getNumInt(io_numIntParms)
     call RuntimeParameters_getNumStr(io_numStrParms)
     call RuntimeParameters_getNumLog(io_numLogParms)

  else
     io_numRealParms = maxParms
     io_numIntParms = maxParms
     io_numStrParms = maxParms
     io_numLogParms = maxParms

     io_numRealScalars = maxParms
     io_numIntScalars = maxParms
     io_numStrScalars = maxParms
     io_numLogScalars = maxParms

     call RuntimeParameters_get('CPNumber', io_checkpointFileNumber)
     call IO_readCheckpoint()
     
  end if

  color = Io_globalMe/io_procGroup
  key = mod(io_globalMe,gr_procGrid(1))
  io_group=color
  call MPI_Comm_split(io_globalComm,color,key,io_comm,ierr)

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  lsize=(blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1)*&
        (blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1)*&
        (blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1)
  start=key*lsize
  gize=lsize*io_procGroup
  arrayDim=1
  call MPI_TYPE_CREATE_SUBARRAY(arrayDim,gsize,lsize,start,&
          MPI_ORDER_FORTRAN,FLASH_REAL,io_datatype,ierr)
  call MPI_TYPE_COMMIT(io_datatype,ierr)

  !Set the initial state of io_outputInStack.  This variable is used by 
  !io_restrictBeforeWrite.F90.  It is .false. for most of the simulation, 
  !and is only briefly .true. when an IO_output routine is in the call stack.
  io_outputInStack = .false.

  return
end subroutine IO_init
