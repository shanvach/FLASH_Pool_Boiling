!!****if* source/Simulation/SimulationMain/unitTest/SolidBody/MultipleWithTriangles/MovingBodies/NoComm/Driver_evolveFlash
!!
!!
!! NAME
!!
!!  Driver_evolveFlash
!!
!! SYNOPSIS
!!
!!  Driver_evolveFlash ()
!!
!! DESCRIPTION
!!
!!
!! NOTES
!!
!!  The Driver unit uses a few unit scope variables that are
!!  accessible to all routines within the unit, but not to the
!!  routines outside the unit. These variables begin with "dr_"
!!  like, dr_globalMe or dr_dt, dr_beginStep, and are stored in fortran
!!  module Driver_data (in file Driver_data.F90. The other variables
!!  are local to the specific routine and do not have the prefix "dr_"
!!
!!
!!***

subroutine Driver_evolveFlash ()

  use Driver_data,         ONLY : dr_globalMe, dr_nbegin,         &
                                  dr_nstep, dr_dt, dr_simTime,    &
                                  dr_nend, dr_dtOld, dr_dtNew,    &
                                  dr_tmax, dr_wallClockTimeLimit, &
                                  dr_elapsedWCTime
  use Driver_interface,    ONLY : Driver_computeDt,               &
                                  Driver_getElapsedWCTime
  use gr_sbInterface,      ONLY : gr_sbCreateParticles
  use gr_sbData,           ONLY : gr_sbBodyInfo, gr_sbNumBodies
  use Logfile_interface,   ONLY : Logfile_stamp,                 &
                                  Logfile_close
  use Timers_interface,    ONLY : Timers_start,                  &
                                  Timers_stop,                   &
                                  Timers_getSummary
  implicit none

# include "constants.h"
# include "Flash.h"
include "Flash_mpi.h"

  integer :: temp,i, ierr, b
  
  character (len=MAX_STRING_LENGTH) :: step

  call Logfile_stamp ('Entering evolution routine' , '[Driver_evolveFlash]')
  call Timers_start  ("evolution")

  call Driver_computeDt(dr_nbegin,  dr_nstep, &
       dr_simTime, dr_dtOld, dr_dtNew)
  
  ! store new                                                                                                                                                                
  dr_dt = dr_dtNew

  !Step forward in time
  do dr_nstep = dr_nbegin, dr_nend
     dr_simTime = dr_simTime + dr_dt

     write(step, '(i10)' ) dr_nstep
     call Logfile_stamp("Evolution step "//step)
     call log_memory_usage()

     call Timers_start("Particles_advance")   
     call gr_ptAdvance(dr_dtOld, dr_dt)
     call Timers_stop("Particles_advance")
     
     call Timers_start("Grid_getBoundboxCentroids")
     call Grid_getBoundboxCentroids()
     call Timers_stop("Grid_getBoundboxCentroids")

     call Grid_sbSelectMaster()
     call gr_sbCreateParticles()
     call Grid_sbBroadcastParticles()
     call gr_sbFinalize()
     
     ! Compute next step dt:
     ! backup needed old
     dr_dtOld = dr_dt
     
     ! calculate new
     call Driver_computeDt(dr_nbegin,  dr_nstep,      &
          dr_simTime, dr_dtOld, dr_dtNew)
     ! store new
     dr_dt = dr_dtNew
        
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
  
  call Timers_stop   ("evolution")
  call Logfile_stamp ('Exiting evolution routine' , '[Driver_evolveFlash]')

  call Timers_getSummary (dr_nstep)

  call Logfile_close ()
  return
end subroutine Driver_evolveFlash


subroutine log_memory_usage()
  use driver_data, ONLY : dr_globalMe, dr_globalComm
  use ut_sysMemInterface, ONLY : ut_sysMemInit, ut_sysMemSummaryStats
  use ut_sysMemData, ONLY : memsummary_t
  use Logfile_interface, ONLY: Logfile_stamp
  implicit none
  type(memsummary_t), allocatable, save, dimension(:) :: memSummary
  integer :: i, numMemMeasurements
  character (len=12) :: minMem_string, maxMem_string, avgMem_string

  call ut_sysMemInit(dr_globalComm, "all", numMemMeasurements)
  allocate(memSummary(numMemMeasurements))

  call ut_sysMemSummaryStats(numMemMeasurements, memSummary)
  do i = 1, numMemMeasurements
     if (dr_globalMe == MASTER_PE) then
        write (minMem_string, '(f12.2)') memSummary(i) % min
        write (maxMem_string, '(f12.2)') memSummary(i) % max
        write (avgMem_string, '(f12.2)') memSummary(i) % avg

        call Logfile_stamp(trim(memSummary(i) % description) // & 
             trim(minMem_string) // ' (min)  ' // &
             trim(maxMem_string) // ' (max)  ' // &
             trim(avgMem_string) // ' (avg) ', tag="memory")
     end if
  end do

  deallocate(memSummary)
end subroutine log_memory_usage
