!!****if* source/Simulation/SimulationMain/unitTest/PFFT_XYperZneu_MC/Driver_evolveFlash
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


#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif

subroutine Driver_evolveFlash()

#include "Flash.h"

  use Driver_data, ONLY: dr_globalMe, dr_globalComm, dr_nbegin,       &
                         dr_nend, dr_dt, dr_wallClockTimeLimit, &
                         dr_tmax, dr_simTime, dr_redshift,      &
                         dr_nstep, dr_dtOld, dr_dtNew,          &
                         dr_restart, dr_elapsedWCTime

  use Driver_interface, ONLY : Driver_sourceTerms, Driver_computeDt, &
       Driver_getElapsedWCTime
  use Logfile_interface,ONLY : Logfile_stamp, Logfile_close
  use Timers_interface, ONLY : Timers_start, Timers_stop, &
                               Timers_getSummary
  use Particles_interface, ONLY : Particles_advance, Particles_dump

  use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                                Grid_getListOfBlocks, &
                                Grid_updateRefinement, &
                                Grid_getBlkPtr, &
                                Grid_getBlkIndexLimits, &
                                Grid_solvePoisson, &
                                Grid_releaseBlkPtr


  use IO_interface,      ONLY : IO_output,IO_outputFinal
   
  use gr_interface, ONLY : gr_findMean
  
  !use Grid_data, ONLY : gr_meshMe


#ifdef FLASH_GRID_PARAMESH
  use tree, ONLY: lrefine, lrefine_max
#endif


  implicit none

#include "constants.h"

#ifdef FLASH_GRID_UG
#else
#include "Multigrid.h"
#endif
 !include "Multigrid.h"
 include "Flash_mpi.h"

  integer   :: localNumBlocks

  integer :: blockCount
  integer :: blockList(MAXBLOCKS)
  integer :: sweepDummy

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, pointer, dimension(:,:,:,:) :: solnData

  ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(3,2) :: strBuff
  character(len=15) :: numToStr

  logical :: endRun  

  real :: turbkin
  real, parameter :: diss = 0.4

  integer, dimension(6) :: bc_types 
  real, dimension(2,6)  :: bc_values
  real poisfact

  integer blkpoints, blkpointsaux
  real L2_err, L2_erraux, Linf_err, Linf_erraux
  
  integer blockID,lb,ierr

  integer TA(2),count_rate
  real   ET
  
  real meanPHI,meanPHIaux,meanANL,meanANLaux 

  real, dimension(2,MDIM) :: boundBox



  integer ii,jj,kk
  real :: coord(MDIM),bsize(MDIM)



  meanPHI = 0.0
  meanPHIaux =0.0 
  meanANL = 0.0 
  meanANLaux = 0.0

  call Logfile_stamp( 'Entering evolution loop' , '[Driver_evolveFlash]')
  call Timers_start("evolution")

  ! Boundary Conditions:
#ifdef FLASH_GRID_UG

  bc_types(1:2) = PERIODIC
  bc_types(3:4) = PERIODIC
  bc_types(5:6) = OUTFLOW  
  bc_values = 0.

#else
  bc_types(1:2) = MG_BND_PERIODIC
  bc_types(3:4) = MG_BND_PERIODIC
  bc_types(5:6) = MG_BND_PERIODIC  
  bc_values = 0.


#endif


  ! Call Poisson Solver
  if (dr_globalMe .eq. 0) write(*,*) 'Into Grid Solve Poisson ..'


  call mpi_barrier(MPI_COMM_WORLD,ierr)
  if (dr_globalMe .eq. 0) CALL SYSTEM_CLOCK(TA(1),count_rate)  


  !!! Write files: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  call Grid_getLocalNumBlks(lnblocks2)
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

!!$  if(gr_meshMe .eq. 0) then
!!$  open(unit=33, file="./IOData/BlockSolMG1_3D.dat", status='replace')
!!$  else
!!$  open(unit=33, file="./IOData/BlockSolMG2_3D.dat", status='replace')
!!$  endif
!!$  write(*,*) 'Writing BlockSolMG_3D.dat ...'
!!$  write(33,*) NGUARD,blockCount
!!$  write(33,*) NXB,NYB,NZB
!!$
!!$  do lb = 1,blockCount
!!$
!!$     blockID = blockList(lb)
!!$     ! Point to blocks center and face vars:
!!$     call Grid_getBlkPtr(blockID,solnData,CENTER)
!!$
!!$     call Grid_getBlkBoundBox(blockId,boundBox)
!!$     bsize(1:NDIM) = boundBox(2,1:NDIM) - boundBox(1,1:NDIM)
!!$
!!$     call Grid_getBlkCenterCoords(blockId,coord)
!!$
!!$     ! Get Blocks internal limits indexes:
!!$     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
!!$
!!$
!!$     write(33,*) coord(1:NDIM)
!!$     write(33,*) bsize(1:NDIM)
!!$
!!$
!!$     do kk = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)    
!!$     do jj = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
!!$     do ii = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
!!$
!!$           write(33,*) solnData(VSRC_VAR,ii,jj,kk) 
!!$
!!$     enddo
!!$     enddo
!!$     enddo
!!$
!!$     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
!!$     
!!$
!!$  enddo
!!$  close(33);  write(*,*) 'BlockSolMG_3D.dat written.'
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1  
!!$  call Driver_abortFlash("Finished writing files")


  poisfact = 1.0 
  call Grid_solvePoisson(VPHI_VAR, VSRC_VAR, bc_types, bc_values, poisfact)


  call mpi_barrier(MPI_COMM_WORLD,ierr)
  if (dr_globalMe .eq. 0) then
      CALL SYSTEM_CLOCK(TA(2),count_rate)
      ET=REAL(TA(2)-TA(1),8)/count_rate
      write(*,*) 'Poisson Solver time =',ET
   endif


  

  call Timers_stop("evolution")
  call Logfile_stamp( 'Exiting evolution loop' , '[Driver_evolveFlash]')


  ! Mean Phi
  call gr_findMean(VPHI_VAR,2,.false.,meanPHI)
  call gr_findMean(VANL_VAR,2,.false.,meanANL)

  !!Step forward in time. See bottom of loop for time step calculation.
  call Grid_getLocalNumBlks(localNumBlocks)
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  ! Check error in the solution:
  L2_err = 0.
  blkpoints = 0
  Linf_err = 0.
!!$  meanPHI = 0.
!!$  meanANL = 0.
  do lb = 1,blockCount

     blockID = blockList(lb)

     if(lrefine(blockID) .eq. lrefine_max) then

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC) 

     blkpoints = blkpoints + &
                 (blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1) * &
                 (blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1) * &
                 (blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1)


!!$     meanPHI = meanPHI + &
!!$               sum( solnData(VPHI_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),     &
!!$                                      blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),     & 
!!$                                      blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))) !/real(blkpoints)
!!$
!!$     meanANL = meanANL + &
!!$               sum( solnData(VANL_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),     &
!!$                                      blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),     & 
!!$                                      blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))) !/real(blkpoints)

     solnData(VERR_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),           &
                       blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),           & 
                       blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))   =       &
     abs(  solnData(VPHI_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),     &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),     & 
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))   - &
           solnData(VANL_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),     &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),     & 
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))   )
     !solnData(VERR_VAR,:,:,:) = solnData(VERR_VAR,:,:,:) - meanPHI
     ! L2 norm of error:
     L2_err = L2_err + sum( solnData(VERR_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),           &
                                              blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),           & 
                                              blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))**2.)

     ! Linf norm of error:
     Linf_err = max(Linf_err,maxval(solnData(VERR_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),   &
                                              blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),           & 
                                              blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) ))
     
     ! Release Pointer
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)

     endif

  enddo  


!!$  ! Sum processors meanPHI
!!$  meanPHIaux = meanPHI
!!$  call MPI_Allreduce(meanPHIaux, meanPHI, 1, FLASH_REAL,&
!!$                     MPI_SUM, MPI_COMM_WORLD, ierr)
!!$
!!$  ! Sum processors meanANL
!!$  meanANLaux = meanANL
!!$  call MPI_Allreduce(meanANLaux, meanANL, 1, FLASH_REAL,&
!!$                     MPI_SUM, MPI_COMM_WORLD, ierr)
  
  ! Sum processors L2 norm of error squared
  blkpointsaux = blkpoints
  call MPI_Allreduce(blkpointsaux, blkpoints, 1, FLASH_INTEGER,&
                     MPI_SUM, MPI_COMM_WORLD, ierr)

  ! Sum processors L2 norm of error squared
  L2_erraux = L2_err
  call MPI_Allreduce(L2_erraux, L2_err, 1, FLASH_REAL,&
                     MPI_SUM, MPI_COMM_WORLD, ierr)

  ! Compute L2 norm for whole domain
  L2_err = sqrt(1./real(blkpoints) * L2_err)

  ! Sum processors Linf norm of error
  Linf_erraux = Linf_err
  call MPI_Allreduce(Linf_erraux, Linf_err, 1, FLASH_REAL,&
                     MPI_MAX, MPI_COMM_WORLD, ierr)

  if ( dr_globalMe .eq. 0) then
     write(*,*) 'Mean Anl, Num=',meanANL,meanPHI     
     write(*,*) 'L2 error = ',L2_err
     write(*,*) 'Linf error = ',Linf_err
     write(*,*) '############################################'

  endif

  ! Export to Tecplot:
  call outtotecplot(dr_globalMe,0.0,1.,1,0,0.0,blockList,blockCount,0)



!  if(.NOT.endRun) call IO_outputFinal( dr_globalComm)
  !call Timers_getSummary(dr_globalMe, dr_nstep)
  call Logfile_stamp( "FLASH run complete.", "LOGFILE_END")
  call Logfile_close()

  return
  
end subroutine Driver_evolveFlash



