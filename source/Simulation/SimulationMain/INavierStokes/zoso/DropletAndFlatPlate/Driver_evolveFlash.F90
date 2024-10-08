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


#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif

#define WRITE_TO_TECPLOT 1

subroutine Driver_evolveFlash()

  use Driver_data, ONLY: dr_globalMe, dr_nbegin,                    &
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
  use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                                Grid_getListOfBlocks, &
                                Grid_updateRefinement

  use IO_interface,      ONLY : IO_output,IO_outputFinal

  use IO_data , ONLY : IO_checkpointFileIntervalStep, io_plotFileNumber, IO_plotFileIntervalTime, IO_plotFileIntervalStep

  use tree, only : grid_changed

  use IncompNS_data, only : ins_cflflg
  use Multiphase_interface, only: Multiphase

  implicit none

#include "constants.h"
#include "Flash.h"
 include "Flash_mpi.h"

  integer   :: localNumBlocks

  integer :: blockCount
  integer :: blockList(MAXBLOCKS)
  integer :: sweepDummy
  
  ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(3,2) :: strBuff
  character(len=15) :: numToStr

  logical :: endRun, gridChanged  


  integer :: count, firstfileflag

  logical :: tecplot_flg

  integer :: ibd,NumBodies
!-----------------------------------------------------------------------------------------

!KPD
if (dr_nstep .eq. 1) grid_changed = 1

!-----------------------------------------------------------------------------------------
  endRun = .false.

  call Logfile_stamp( 'Entering evolution loop' , '[Driver_evolveFlash]')
  call Timers_start("evolution")

  ! Initial Timestep:
  ! backup needed old
  dr_dtOld = dr_dt

  ! calculate new
  call Driver_computeDt( dr_nbegin,  dr_nstep,      &
                        dr_simTime, dr_dtOld, dr_dtNew)
  ! store new
  dr_dt = dr_dtNew

  !print *,"dr_dt: ",dr_dt

  !write(*,*) 'dr_dt ===',dr_dt

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  !- kpd - Restart Plot File #
  !---------------------------
  if (dr_restart .eqv. .TRUE.) then
     !count = io_plotFileNumber
     count = io_plotFileNumber -1 ! Shizhao
  else
     count = 0
  end if

  firstfileflag = 0
  call outtotecplot(dr_globalMe,dr_simtime,dr_dt,dr_nstep,count, &
                    0.0,blockList,blockCount,firstfileflag)

  call outtotecplot_uv(dr_globalME,dr_simtime,dr_dt,dr_nstep,count, &
                       0.0,blockList,blockCount,firstfileflag)

  firstfileflag = 1
  do dr_nstep = dr_nBegin, dr_nend
     
     !!Step forward in time. See bottom of loop for time step calculation.
     call Grid_getLocalNumBlks(localNumBlocks)
     call Grid_getListOfBlocks(LEAF,blockList,blockCount)

     if (dr_globalMe == MASTER_PE) then

        write (numToStr(1:), '(I10)') dr_nstep
        write (strBuff(1,1), "(A)") "n"
        write (strBuff(1,2), "(A)") trim(adjustl(numToStr))
        
        write (numToStr(1:), "(1PE12.6)") dr_simTime
        write (strBuff(2,1), "(A)") "t"
        write (strBuff(2,2), "(A)") trim(adjustl(numToStr))
        
        write (numToStr(1:), "(1PE12.6)") dr_dt
        write (strBuff(3,1), "(A)") "dt"
        write (strBuff(3,2), "(A)") trim(adjustl(NumToStr))
        
        !call Logfile_stamp( strBuff, 3, 2, "step")
     end if


     !--------------------------------------------------------------------
     !--------------------------------------------------------------------
     !- Start Physics Sequence
     !--------------------------------------------------------------------
     !--------------------------------------------------------------------
     !----
     dr_simTime = dr_simTime + dr_dt

     call Timers_start("Multiphase")
     call Multiphase( blockCount, blockList,   &
              dr_simTime, dr_dt, dr_dtOld,  sweepDummy,1)
     call Timers_stop("Multiphase")

     call Timers_start("Multiphase")
     call Multiphase( blockCount, blockList,   &
              dr_simTime, dr_dt, dr_dtOld,  sweepDummy,0)
     call Timers_stop("Multiphase")

#ifdef DEBUG_DRIVER
     print*, 'going into IncompNS'
#endif
     !print *,"Time step: ",dr_dt

     call Timers_start("IncompNS")
     call IncompNS( blockCount, blockList,   &
              dr_simTime, dr_dt, dr_dtOld,  sweepDummy)
     call Timers_stop("IncompNS")

#ifdef DEBUG_DRIVER
  print*, 'return from IncompNS timestep'
#endif

     if (dr_globalMe .eq. MASTER_PE) then
        write(*,*) ' '        
        write(*,'(I6,A,g16.8,A,g16.8)') dr_nstep,&
                ', TimeStep= ',dr_dt,', SimTime= ', dr_simTime
     endif     

     if (dr_globalMe .eq. MASTER_PE) &
     write(*,*) '###############################################################################'

     if (ins_cflflg .eq. 1) then ! Constant cfl       
       if (dr_nstep .gt. 1) then
       tecplot_flg = (1/IO_plotFileIntervalTime*MOD(dr_simtime,IO_plotFileIntervalTime) .le. &
                      dr_dt/IO_plotFileIntervalTime)
       else
       tecplot_flg = .false.
       endif
     else                        ! Constant timestep
       tecplot_flg = (MOD(dr_nstep,IO_plotFileIntervalStep) .eq. 0)
     endif

     if (tecplot_flg) then
        ! Write to Grid to Tecplot:
        count = count + 1
        call outtotecplot(dr_globalMe,dr_simtime,dr_dt,dr_nstep,count, &
                          0.0,blockList,blockCount,firstfileflag)

        call outtotecplot_uv(dr_globalME,dr_simtime,dr_dt,dr_nstep,count, &
                       0.0,blockList,blockCount,firstfileflag)

     if (count .gt. 0) firstfileflag = 1
     endif
     !--------------------------------------------------------------------     

     !----
     !--------------------------------------------------------------------
     !--------------------------------------------------------------------
     !- End Physics Sequence
     !--------------------------------------------------------------------
     !--------------------------------------------------------------------

     call Grid_updateRefinement(dr_nstep, dr_simTime, gridChanged )

     ! backup needed old
     dr_dtOld = dr_dt

     ! calculate new
     call Driver_computeDt( dr_nbegin,  dr_nstep,      &
                           dr_simTime, dr_dtOld, dr_dtNew)
     ! store new
     dr_dt = dr_dtNew
    
     call Timers_start("io")
!     call IO_output(dr_simTime,dr_dt,dr_nstep+1,dr_nbegin,endRun)
     call Timers_stop("io")

     if(endRun) exit

     !! the simulation ends before nend iterations if
     !!  (i)   the simulation time is greater than the maximum time (tmax)
     !!  (ii)  the redshift falls below the minimum redshift  
     !!        (also called zfinal)
     !!  (iii) the wall clock time is greater than the maximum 
     !!        (wall_clock_time_max)

     if (dr_simTime >= dr_tmax) then
        if(dr_globalMe == MASTER_PE) then
           print *, "exiting: reached max SimTime"
        endif
        exit
     end if
     
     call Driver_getElapsedWCTime(dr_elapsedWCTime)
     if (dr_elapsedWCTime >  dr_wallClockTimeLimit) then
        if(dr_globalMe == MASTER_PE) then
           print *, "exiting: reached max wall clock time"
        endif
        exit
     end if

  enddo

  call Timers_stop("evolution")
!  if(.NOT.endRun) call IO_outputFinal()
  call Timers_getSummary(dr_nstep)

  return
  
end subroutine Driver_evolveFlash
