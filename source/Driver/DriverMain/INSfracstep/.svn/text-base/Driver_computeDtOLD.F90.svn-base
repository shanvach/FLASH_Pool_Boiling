!!****if* source/Driver/DriverMain/INSfracstep/Driver_computeDt
!!
!! NAME
!!
!!  Driver_computeDt
!!
!!
!! SYNOPSIS
!!
!!  Driver_computeDt( integer(IN) :: nbegin,
!!                    integer(IN) :: nstep, 
!!                    real(IN)    :: simTime
!!                    real(IN)    :: dtOld,
!!                    real(OUT)   :: dtNew )
!!
!! DESCRIPTION
!!
!!  Determine the stability-limited time step.
!!  This timestep is determined using information from the included
!!  physics modules - many different timestep limiters are polled.
!!
!!  The global driver might use a different (hopefully smaller) time
!!  step, to match a file write time (tplot or trstr) or if the
!!  simulation end time has been reached; such possibilities are
!!  not considered here.
!!
!! ARGUMENTS
!!  nbegin   - first step of the simulation (beginStep is only used
!!             to determine if a label header should be written to
!!             the screen)
!!  nstep    - current step of the simulation
!!  simTime  - current simulation time of the run
!!  dtOld    - the dt from the timestep that we just finished 
!!             (it's old because we be using dtOld to calculate 
!!             and return the dt for the next timestep (dtNew)
!!  dtNew    - returned value of the dt calculated for the next timestep
!!
!!
!!
!! NOTES 
!!
!!  variables that begin with "dr_" like, dr_dtMin or dr_dtMax, 
!!  are stored in the data fortran module for the Driver unit, 
!!  Driver_data.  The "dr_" is meant to indicate that the 
!!  variable belongs to the Driver Unit.
!!  all other normally named variables i, j, etc are local variables.
!!
!!
!!*** 

#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif

subroutine Driver_computeDt( nbegin,    nstep,     &
                             simTime,   dtOld, dtNew)


  use Grid_interface, ONLY : Grid_getListOfBlocks,  &
                             Grid_getBlkIndexLimits,&
                             Grid_getDeltas,        &
                             Grid_getBlkPtr,        &
                             Grid_releaseBlkPtr

  use IncompNS_data, ONLY : ins_cflflg, ins_dtspec

  use IncompNS_interface, ONLY : IncompNS_computeDt
  use Driver_data, ONLY : dr_globalMe

  implicit none

#include "constants.h"
#include "Flash.h"
#include "IncompNS.h"
 include "Flash_mpi.h"

  integer, intent(IN) :: nbegin, nstep
  real,    intent(IN) :: simTime    !! current simulation time
  real,    intent(IN) :: dtOld      !! last time step we used
  real,    intent(OUT):: dtNew      !! the new timestep we get. to be returned.
 

  ! Local variables and functions
  integer :: i, j, error, blockID, numLeafBlocks, iout, istat

  real, PARAMETER :: MAX_TSTEP = huge(1.0)

  
  real    :: dtModule, dtLocal

  integer, dimension(MAXBLOCKS) :: blockList

  !!prepatory data structures for passing coords to timestep routines
  real, dimension(MDIM) :: del
  integer, dimension(MDIM) :: index


  !!arrays which hold the starting and ending indices of a block
  integer,dimension(2,MDIM)::blkLimits,blkLimitsGC

  !!coordinate infomration to be passed into physics  
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData
  integer :: isize,jsize,ksize

  integer :: ii,jj
  
  !!Initialize all timestep variables.
  dtLocal = MAX_TSTEP
  
  !! Loop over all local leaf-node blocks
  call Grid_getListOfBlocks(LEAF,blockList, numLeafBlocks)

  do i = 1, numLeafBlocks

     !!Get the coordinate information for all the
     call Grid_getBlkIndexLimits(blockList(i),blkLimits,blkLimitsGC)
     isize = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
     jsize = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
     ksize = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

     blockID=blockList(i)

     call Grid_getDeltas(blockID, del)

     call Grid_getBlkPtr(blockID, facexData,FACEX)
     call Grid_getBlkPtr(blockID, faceyData,FACEY)
#if NDIM == 3
     call Grid_getBlkPtr(blockID, facezData,FACEZ)
#endif

#ifdef DEBUG_DRIVER
     print*,'going to call INS timestep'
#endif

     call IncompNS_computeDt (blockID,         & 
                         isize, jsize, ksize,  &
              del(DIR_X),del(DIR_Y),del(DIR_Z),&
                         blkLimits,blkLimitsGC,&
                         facexData,faceyData,  &
                         facezData,            &
                         dtLocal )

     call Grid_releaseBlkPtr(blockID, facexData, FACEX)
     call Grid_releaseBlkPtr(blockID, faceyData, FACEY)
#if NDIM == 3
     call Grid_releaseBlkPtr(blockID, facezData, FACEZ) 
#endif

#ifdef DEBUG_DRIVER
     print*,'returned from INS timestep'
#endif

  enddo

  ! Find the minimum timestep across all processors and all modules.
  call MPI_AllReduce (dtLocal, dtModule, 1, & 
       MPI_DOUBLE_PRECISION, MPI_MIN, MPI_Comm_World, error)

  dtNew = HUGE(1.0)          ! dt will hold the minimum timestep  

  if (dtModule < dtNew) then
        dtNew = dtModule
  endif

  ! If running at constant timestep check that is smaller that cfl dt: 
  if (ins_cflflg .eq. 0) then
     if( dtNew .ge. ins_dtspec) then
     dtNew = ins_dtspec
     else
        if (dr_globalMe == MASTER_PE) then
           write(*,'(A,g16.8,A,g16.8)') 'SPECIFIED DT= ',ins_dtspec,', LARGER THAN cfl DT= ',dtNew
           write(*,'(A)') 'USING cfl DT ..'
        endif
     endif
  endif
   
  return

end subroutine Driver_computeDt
