!!****if* source/Simulation/SimulationMain/unitTest/BiCG_XYdir_2D/Driver_evolveFlash
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

#include "Flash.h"

subroutine Driver_evolveFlash()

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

  use Grid_interface,    ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
                                GRID_PDE_BND_DIRICHLET, Grid_getLocalNumBlks, &
                                Grid_getListOfBlocks, &
                                Grid_updateRefinement, &
                                Grid_getBlkPtr, &
                                Grid_getBlkIndexLimits, &
                                Grid_solvePoisson, &
                                Grid_releaseBlkPtr

  use gr_interface, ONLY : gr_findMean

  use IO_interface,      ONLY : IO_output,IO_outputFinal

  use RuntimeParameters_interface, ONLY: RuntimeParameters_get, &
       RuntimeParameters_mapStrToInt

   
#ifdef FLASH_GRID_PARAMESH  
 use tree, only : lrefine,lrefine_max
! use gr_mgPfftData, ONLY : gr_mgbcTypes
#endif

  
  implicit none

#include "constants.h"
#ifdef FLASH_GRID_UG
#else
!#include "Multigrid.h"
#endif
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

  integer, dimension(6) :: bc_types, gr_PoissonBcTypes 
  real, dimension(2,6)  :: bc_values
  real poisfact

  integer blkpoints, blkpointsaux
  real L2_err, L2_erraux, Linf_err, Linf_erraux
  
  integer blockID,lb,ierr

  integer TA(2),count_rate
  real   ET
  
  real meanPHI,meanPHIaux,meanANL,meanANLaux 

  integer :: eachBoundary
  character(len=MAX_STRING_LENGTH), dimension(2*MDIM) :: bcTypeStr

  meanPHI = 0.0
  meanPHIaux =0.0 
  meanANL = 0.0 
  meanANLaux = 0.0

  call Logfile_stamp( 'Entering evolution loop' , '[Driver_evolveFlash]')
  call Timers_start("evolution")

  !Initialise to -1 to help us spot errors.
  bc_types = -1

  ! Retrieve Poisson solution Boundary Conditions for each face: 
  call RuntimeParameters_get("xl_boundary_type", bcTypeStr(1))
  call RuntimeParameters_get("xr_boundary_type", bcTypeStr(2))
  if (NDIM >= 2) then
     call RuntimeParameters_get("yl_boundary_type", bcTypeStr(3))
     call RuntimeParameters_get("yr_boundary_type", bcTypeStr(4))
  endif
  if (NDIM == 3) then
     call RuntimeParameters_get("zl_boundary_type", bcTypeStr(5))
     call RuntimeParameters_get("zr_boundary_type", bcTypeStr(6))
  endif


  ! Boundary Conditions:
#ifdef FLASH_GRID_UG

  do eachBoundary = 1, 2*NDIM

     call RuntimeParameters_mapStrToInt(bcTypeStr(eachBoundary), bc_types(eachBoundary))

  enddo

!  bc_types(1:2) = DIRICHLET 
!  bc_types(3:4) = DIRICHLET
!  bc_types(5:6) = DIRICHLET
  bc_values = 0.

#else

  do eachBoundary = 1, 2*NDIM

     call RuntimeParameters_mapStrToInt(bcTypeStr(eachBoundary), gr_PoissonBcTypes(eachBoundary))

     select case(gr_PoissonBcTypes(eachBoundary))
     case (PERIODIC) 
        bc_types(eachBoundary)  = GRID_PDE_BND_PERIODIC
     case (OUTFLOW)
        bc_types(eachBoundary)  = GRID_PDE_BND_NEUMANN
     case (DIRICHLET)
        bc_types(eachBoundary)  = GRID_PDE_BND_DIRICHLET
     case default
           write(*,*) 'In DriverEvolveFlash: Poisson Problem Boundary Condition not supported'
           call Driver_abortFlash('BCs unsupported')        
     end select
  enddo


!  bc_types(1:6) = gr_mgbcTypes(1:6)
!  bc_types(1:2) = MG_BND_DIRICHLET
!  bc_types(3:4) = MG_BND_DIRICHLET
!  bc_types(5:6) = MG_BND_DIRICHLET
  bc_values = 0.


#endif


  ! Call Poisson Solver
  if (dr_globalMe .eq. 0) write(*,*) 'Into Grid Solve Poisson ..'


  call mpi_barrier(dr_globalComm,ierr)
  if (dr_globalMe .eq. 0) CALL SYSTEM_CLOCK(TA(1),count_rate)  


  poisfact = 1.0 
  call Grid_solvePoisson(VPHI_VAR, VSRC_VAR, bc_types, bc_values, poisfact)


  call mpi_barrier(dr_globalComm,ierr)
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
  do lb = 1,blockCount

     blockID = blockList(lb)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC) 


     solnData(VERR_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),           &
                       blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),           & 
                       blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))   =       &
     abs(  solnData(VPHI_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),     &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),     & 
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))   - &
           solnData(VANL_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),     &
                             blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),     & 
                             blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))   )

#ifdef FLASH_GRID_PARAMESH
     if (lrefine(blockID) .eq. lrefine_max) then
#endif

     blkpoints = blkpoints + &
                 (blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1) * &
                 (blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1) * &
                 (blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1)

     ! L2 norm of error:
     L2_err = L2_err + sum( solnData(VERR_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),           &
                                              blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),           & 
                                              blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))**2.)

     ! Linf norm of error:
     Linf_err = max(Linf_err,maxval(solnData(VERR_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),   &
                                              blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),           & 
                                              blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) ))
#ifdef FLASH_GRID_PARAMESH
     endif
#endif

     ! Release Pointer
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)


  enddo  
  
  ! Sum processors L2 norm of error squared
  blkpointsaux = blkpoints
  call MPI_Allreduce(blkpointsaux, blkpoints, 1, FLASH_INTEGER,&
                     MPI_SUM, dr_globalComm, ierr)

  ! Sum processors L2 norm of error squared
  L2_erraux = L2_err
  call MPI_Allreduce(L2_erraux, L2_err, 1, FLASH_REAL,&
                     MPI_SUM, dr_globalComm, ierr)

  ! Compute L2 norm for whole domain
  L2_err = sqrt(1./real(blkpoints) * L2_err)

  ! Sum processors Linf norm of error
  Linf_erraux = Linf_err
  call MPI_Allreduce(Linf_erraux, Linf_err, 1, FLASH_REAL,&
                     MPI_MAX, dr_globalComm, ierr)

  if ( dr_globalMe .eq. 0) then
     write(*,*) 'Mean Anl, Num=',meanANL,meanPHI     
     write(*,*) 'L2 error = ',L2_err
     write(*,*) 'Linf error = ',Linf_err
     write(*,*) '############################################'

  endif

  ! Export to Tecplot:
  call outtotecplot(dr_globalMe,0.0,1.,1,0,0.0,blockList,blockCount,0)

!  if(.NOT.endRun) call IO_outputFinal()
  call Timers_getSummary( dr_nstep)
  call Logfile_stamp( "FLASH run complete.", "LOGFILE_END")
  call Logfile_close()

  return
  
end subroutine Driver_evolveFlash



