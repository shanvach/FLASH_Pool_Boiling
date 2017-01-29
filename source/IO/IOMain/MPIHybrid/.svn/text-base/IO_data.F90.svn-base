!!****if* source/IO/IOMain/MPIHybrid/IO_data
!!
!! NAME
!!  IO_data
!!
!! SYNOPSIS
!!
!!  IO_data()
!!
!! DESCRIPTION 
!!  
!!  Holds all the IO data that is needed by the IO unit     
!!
!! ARGUMENTS
!!
!!  none    
!!
!!
!!***

Module IO_data
  integer,save :: io_comm,io_me,io_datatype,io_procGroup,io_group
  integer, save :: io_globalMe, io_globalNumProcs
  integer, save :: io_ptNumber, io_checkpointFileNumber 
  integer, save :: io_nextStep !DEV: verify we are using this
  integer, save :: io_checkpointFileIntervalStep
  real, save :: io_checkpointFileIntervalTime, io_tplot, io_ptplot, io_zrstrt, io_pzplot 
  real, save :: io_wallClockCheckpoint, io_zpnext, io_zplot
  real, save :: io_tpnext, io_tinitial, io_zinitial, io_trnext, io_zrnext

  integer, save :: io_rollingCheckpoint
  integer, save :: io_memoryStatFreq, io_integralFreq
  character(len=MAX_STRING_LENGTH), save :: io_geometry
  real, save :: io_CPUSeconds, io_lastCPUSeconds, io_lastWallClockCheckpoint
  logical, save :: io_restart, io_bytePack
  integer, parameter :: maxParms = 200

  !these are for the uniform grid
  integer, save :: io_iprocs, io_jprocs, io_kprocs


  character (len=MAX_STRING_LENGTH) :: io_parmStrName = "runtime parameters"
  character (len=MAX_STRING_LENGTH) :: io_scalarStrName = "scalars"


  type (context_type), save :: io_scalar  


  !!These are for the checkpoint file
  character (len=MAX_STRING_LENGTH), save :: io_baseName, io_outputDir

  !!Store the integral quantities file name (the .dat file)
  character (len=MAX_STRING_LENGTH), save :: io_statsFileName


  !!These are for the linked list output
  integer :: io_numRealParms, io_numIntParms, io_numStrParms, io_numLogParms
  integer :: io_numRealScalars, io_numIntScalars, io_numStrScalars, io_numLogScalars

  character (len=MAX_STRING_LENGTH), allocatable, save :: io_intParmNames(:), io_intScalarNames(:)
  integer, allocatable, save :: io_intParmValues(:), io_intScalarValues(:)

  character(len=MAX_STRING_LENGTH),allocatable,save :: io_realParmNames(:), io_realScalarNames(:)
  real, allocatable, save :: io_realParmValues(:), io_realScalarValues(:) 

  character (len=MAX_STRING_LENGTH), allocatable, save :: io_strParmNames(:), io_strScalarNames(:)
  character (len=MAX_STRING_LENGTH), allocatable, save :: io_strParmValues(:), io_strScalarValues(:)

  character (len=MAX_STRING_LENGTH), allocatable, save :: io_logParmNames(:), io_logScalarNames(:)
  logical, allocatable, save :: io_logParmValues(:), io_logScalarValues(:)
  integer, allocatable, save :: io_logToIntParmValues(:), io_logToIntScalarValues(:)

  logical, save :: io_runningParticles = .false.

  integer, save :: io_plotVar(UNK_VARS_BEGIN:UNK_VARS_END), io_nPlotVars

  ! create a temporary array to hold the 4 character variable names
  character (len=4), save :: io_unklabels(UNK_VARS_BEGIN:UNK_VARS_END)
  character (len=4), save :: io_plotVarStr(UNK_VARS_BEGIN:UNK_VARS_END)


  
  character (len=MAX_STRING_LENGTH), save                  :: io_flashRelease = ' '
  character (len=MAX_STRING_LENGTH), save  :: io_buildDate, io_buildDir, io_buildMachine, io_fileCreationTime
  character (len=MAX_STRING_LENGTH), save  :: io_setupTimeStamp, io_buildTimeStamp
  character (len=400), save                :: io_setupCall, io_cflags, io_fflags

  !single precision or double precision

  logical, save :: io_doublePrecision



#ifdef FIXEDBLOCKSIZE

  logical, save :: io_fixedBlockSize = .TRUE.

  integer, save :: io_ilo = GRID_ILO
  integer, save :: io_ihi = GRID_IHI
  integer, save :: io_jlo = GRID_JLO
  integer, save :: io_jhi = GRID_JHI
  integer, save :: io_klo = GRID_KLO
  integer, save :: io_khi = GRID_KHI

  real, allocatable, save :: io_unkBuf(:,:,:,:,:)

#else
  logical, save :: io_fixedBlockSize = .FALSE.
  integer, save :: io_iguard, io_jguard, io_kguard
  integer, save :: io_iloGC, io_ihiGC, io_jloGC, io_jhiGC, iio_kloGC, io_khiGC
  integer, save :: io_ilo, io_ihi, io_jlo, io_jhi, io_klo, io_khi

#endif 

  !Used in data restriction decision.
  logical, save :: io_outputInStack

end module IO_data

