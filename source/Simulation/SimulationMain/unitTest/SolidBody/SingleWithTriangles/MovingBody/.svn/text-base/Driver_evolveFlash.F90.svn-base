!!****if* source/Simulation/SimulationMain/unitTest/SolidBody/SingleWithTriangles/MovingBody/Driver_evolveFlash
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
  use gr_sbInterface,      ONLY : gr_sbCreateGroups, gr_sbCreateParticles
  use gr_sbData,           ONLY : gr_sbBodyInfo
  use Logfile_interface,   ONLY : Logfile_stamp,                 &
                                  Logfile_close
  use Timers_interface,    ONLY : Timers_start,                  &
                                  Timers_stop,                   &
                                  Timers_getSummary
  implicit none

# include "constants.h"
# include "Flash.h"
include "Flash_mpi.h"

  character (len=20) :: fileName

  logical :: perfect
  integer :: temp,i, ierr
  
  integer, parameter     :: fileUnit = 2
  integer, dimension (4) :: prNum
!
!
!   ...Give a unique unitTest filename to each processor.
!
!  

  temp = dr_globalMe
  do i = 1,4
     prNum(i)= mod(temp,10)
     temp = temp/10
  end do
  filename = "unitTest_"//char(48+prNum(4))//char(48+prNum(3))//&
                                 char(48+prNum(2))//char(48+prNum(1))
  open  (fileUnit,file=fileName)
  write (fileUnit,'("P",I0)') dr_globalMe
!
!
!   ...Call the opacity unit test routine.
!
!  
  call Logfile_stamp ('Entering evolution routine' , '[Driver_evolveFlash]')
  call Timers_start  ("evolution")

  dr_dtOld = dr_dt

  !calculate new
  call Driver_computeDt(dr_nbegin,  dr_nstep, &
       dr_simTime, dr_dtOld, dr_dtNew)

  ! store new                                                                                                                                                                
  dr_dt = dr_dtNew

  !Step forward in time
  do dr_nstep = dr_nbegin, dr_nend
     dr_simTime = dr_simTime + dr_dt
     call Timers_start("Particles_advance")
     call gr_ptAdvance(dr_dtOld, dr_dt)
     call Timers_stop("Particles_advance")

     call Timers_start("Grid_getBoundboxCentroids")
     call Grid_getBoundboxCentroids()
     call Timers_stop("Grid_getBoundboxCentroids")
     
     call Grid_sbCreateGroups()
     call Grid_sbSelectMaster()
     call gr_sbCreateParticles()
     call Grid_sbBroadcastParticles()
     
     call Grid_solidBodyUnitTest (fileUnit,perfect)

     if (gr_sbBodyInfo(1) % comm /= MPI_COMM_NULL) then
        call MPI_Comm_free(gr_sbBodyInfo(1) % comm, ierr)
     endif

     if (perfect) then
        write (fileUnit,'("All results conformed with expected values.")')
     else
        write (fileUnit,'("Failure in Solid Body unit test")')
     end if
     close (fileUnit)

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
