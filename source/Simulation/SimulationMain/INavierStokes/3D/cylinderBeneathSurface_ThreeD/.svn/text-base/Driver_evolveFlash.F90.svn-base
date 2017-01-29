!!****if* source/Simulation/SimulationMain/INavierStokes/3D/Snorkel_mcHYPRE_VD_wFS/Driver_evolveFlash
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
#define DEDUG_ALL

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
                             Grid_conserveFluxes,    &
                             Grid_conserveField,     &
                             Grid_solvePoisson, Grid_getBlkBoundBox, Grid_getBlkCenterCoords

  use gr_interface, ONLY : gr_findMean

!  use MHD_interface,     ONLY : MHD
!  use Gravity_interface, ONLY :  Gravity_potentialListOfBlocks
  use IO_interface,      ONLY : IO_output,IO_outputFinal

  use IO_data , ONLY : io_plotFileNumber, IO_plotFileIntervalStep, IO_plotFileIntervalTime

!  use Cosmology_interface, ONLY:  Cosmology

  !use ins_interface, only  :  ins_velomg2center, &
  use ins_interface, only  :  ins_divergence, &
                              ins_fluxfix,    &
                              ins_fluxfix_p,  &
                              ins_corrector

  use IncompNS_data, only : ins_cflflg

  use gr_sbData, only : gr_sbNumBodies, gr_sbBodyInfo

  use Grid_data, only : gr_meshMe,gr_meshComm

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

  logical :: endRun

  logical :: tecplot_flg
  integer count, firstfileflag
  
  character(len=6) :: index_ibd,index_count  
  real :: xpt, ypt
  integer :: i,ibd

  integer :: ierr

  logical :: gridChanged

  if (dr_nstep .eq. 1) gridChanged = .TRUE.

  endRun = .false.

  call Logfile_stamp( 'Entering evolution loop' , '[Driver_evolveFlash]')
  call Timers_start("evolution")

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

 ! FIX FLUXES FOR USTAR: (Only for AMR grids)
  ! --- ------ --- -----
!!$#ifdef FLASH_GRID_PARAMESH
!!$
!!$  call Grid_conserveFluxes(dr_MyPE, ALLDIR, -1)
!!$  ! Fix fluxes at block boundaries
!!$  call ins_fluxfix(dr_myPE,NGUARD,nxc,nyc,nzc,nxc-1,nyc-1,nzc-1,&
!!$                   blockCount,blockList)
!!$#endif
!!$  write(*,*) 'AFTER INITIAL FLUX CONS'

  ! Initial Timestep:
  ! backup needed old
  dr_dtOld = dr_dt

  ! calculate new
  call Driver_computeDt(dr_nbegin,  dr_nstep,      &
                        dr_simTime, dr_dtOld, dr_dtNew)
  ! store new
  dr_dt = dr_dtNew

  !write(*,*) 'dr_dt ===',dr_dt

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  count = io_plotFileNumber-1
  firstfileflag = 0 ! io_plotFileNumber
!  dr_nstep = 0
  call outtotecplot(dr_globalMe,dr_simtime,dr_dt,dr_nstep,count, &
                    0.0,blockList,blockCount,firstfileflag)
!  call outtotecplot_uv(dr_globalMe,dr_simtime,dr_dt,dr_nstep,count, &
!                     0.0,blockList,blockCount,firstfileflag)

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

     if (gr_meshMe .eq. MASTER_PE) then
     open(unit=113,file='./IOData/geo'//'9999'//'.plt',form='formatted')
!     write(113,'(A,G12.5,A)')'TEXT X=75,Y=5,F=HELV-BOLD,C=RED,T=" T =',dr_simTime,'"'
     close(113)
     endif

     call MPI_BARRIER(gr_meshComm,ierr)

     do ibd=1,gr_sbNumBodies

        if (gr_sbBodyInfo(ibd)%myPE .eq. gr_sbBodyInfo(ibd)%BodyMaster)then


        open(unit=113,file='./IOData/geo'//'9999'//'.plt',form='formatted', &
             status='old',position='append')

!        call int2char(ibd,index_ibd)
        write(113,'(A)') 'VARIABLES = "X" , "Y", "Z"'
        write(113,'(2(A,I8),A)')                 &
        'ZONE T=Body'//'9999'//', N=',              &
        gr_sbBodyInfo(ibd)%NumVertices,', E=',gr_sbBodyInfo(ibd)%NumAelem,        &
        ',DATAPACKING = POINT, ZONETYPE = FETRIANGLE'

        do i = 1,gr_sbBodyInfo(ibd)%NumVertices
           write(113,'(2E16.8)') gr_sbBodyInfo(ibd)%xb(i),gr_sbBodyInfo(ibd)%yb(i),gr_sbBodyInfo(ibd)%zb(i)
        enddo

        write(113,'(A)') ' '
        Do i=1,gr_sbBodyInfo(ibd)%NumAelem
           write(113,'(3(I8))')                        &
                  gr_sbBodyInfo(ibd)%AELEM(2,i),            &
                  gr_sbBodyInfo(ibd)%AELEM(3,i),gr_sbBodyInfo(ibd)%AELEM(4,i)
        enddo

        close(113)

        end if

     enddo

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

  firstfileflag = 1
!print*,"Calling gr_sbFinalize from Driver_evolveFlash"
!  call gr_sbFinalize() !Added PM
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
        
        call Logfile_stamp( strBuff, 3, 2, "step")
     end if


     !--------------------------------------------------------------------
     !- Start Physics Sequence
     !----
#ifdef DEBUG_DRIVER
     print*, 'going into IncompNS'
#endif
     dr_simTime = dr_simTime + dr_dt

     call Timers_start("IncompNS")
     call IncompNS(blockCount, blockList,   &
                   dr_simTime, dr_dt, dr_dtOld,  sweepDummy)
     call Timers_stop("IncompNS")

#ifdef DEBUG_DRIVER
  print*, 'return from IncompNS timestep'
#endif

!!$     call Timers_start("sourceTerms")
!!$     call Driver_sourceTerms(blockCount, blockList, dr_dt)
!!$     call Timers_stop("sourceTerms")
!!$#ifdef DEBUG_DRIVER
!!$     print*,'done source terms'
!!$     print*, 'return from Drivers_sourceTerms '
!!$#endif
     call Timers_start("Particles_advance")
     call Particles_advance(dr_dtOld, dr_dt)
#ifdef DEBUG_DRIVER
     print*, 'return from Particles_advance '
#endif
     call Timers_stop("Particles_advance")     
!!$     call Gravity_potentialListOfBlocks(blockCount,blockList)
!!$#ifdef DEBUG_DRIVER
!!$     print*, 'return from Gravity_potential '
!!$#endif

     !----
     !- End Physics Sequence
     !--------------------------------------------------------------------

!!$     call Timers_start("Grid_updateRefinement")
     call Grid_updateRefinement( dr_nstep, dr_simTime, gridChanged)
!!$     call Timers_stop("Grid_updateRefinement")

     !----
     !- Output results and plot files: 
     !--------------------------------------------------------------------

     ! Output to Tecplot
     !if (ins_cflflg .eq. 1) then ! Constant cfl       
     if (ins_cflflg .eq. 0) then ! Constant cfl      KPD 
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
     count = count + 1


     call outtotecplot(dr_globalMe,dr_simtime,dr_dt,dr_nstep,count, &
                       0.0,blockList,blockCount,firstfileflag)
!     call outtotecplot_uv(dr_globalMe,dr_simtime,dr_dt,dr_nstep,count, &
!                       0.0,blockList,blockCount,firstfileflag)


     if (gr_meshMe .eq. MASTER_PE) then 
!     call int2char(count,index_count)
     open(unit=113,file='./IOData/geo'//'1234'//'.plt',form='formatted')
     write(113,'(A,G12.5,A)')'TEXT X=75,Y=5,F=HELV-BOLD,C=RED,T=" T =',dr_simTime,'"'
     close(113)
     endif

     call MPI_BARRIER(gr_meshComm,ierr) 

     do ibd=1,gr_sbNumBodies
 
        if (gr_sbBodyInfo(ibd)%myPE .eq. gr_sbBodyInfo(ibd)%BodyMaster)then                
 
        open(unit=113,file='./IOData/geo'//'1234'//'.plt',form='formatted', &
             status='old',position='append')

!        call int2char(ibd,index_ibd)
        write(113,'(A)') 'VARIABLES = "X" , "Y", "Z"'
        write(113,'(2(A,I8),A)')                 &
        'ZONE T=Body'//index_ibd//', N=',              &
        gr_sbBodyInfo(ibd)%NumVertices,', E=',gr_sbBodyInfo(ibd)%NumAelem,        &
        ',DATAPACKING = POINT, ZONETYPE = FETRIANGLE'

        do i = 1,gr_sbBodyInfo(ibd)%NumVertices
           write(113,'(2E16.8)') gr_sbBodyInfo(ibd)%xb(i),gr_sbBodyInfo(ibd)%yb(i),gr_sbBodyInfo(ibd)%zb(i)
        enddo

        write(113,'(A)') ' '
        Do i=1,gr_sbBodyInfo(ibd)%NumAelem
           write(113,'(3(I8))')                        &
                  gr_sbBodyInfo(ibd)%AELEM(2,i),            &
                  gr_sbBodyInfo(ibd)%AELEM(3,i),gr_sbBodyInfo(ibd)%AELEM(4,i)
        enddo

        close(113)

        endif

        call MPI_BARRIER(gr_meshComm,ierr)

     enddo
    

     if (count .gt. 0) firstfileflag = 1
     endif


     !! Average Velocities and Vorticity to cell-centers
     !call ins_velomg2center(blocklist,blockcount) 

     ! Flash Output Routine:
     call Timers_start("io")
     call IO_output(dr_simTime,dr_dt,dr_nstep,dr_nbegin,endRun)
     call Timers_stop("io")

     !--------------------------------------------------------------------


     if (dr_globalMe .eq. MASTER_PE) then
        write(*,*) ' '        
        write(*,'(I6,A,g16.8,A,g16.8)') dr_nstep,&
                ', TimeStep= ',dr_dt,', SimTime= ', dr_simTime
     endif     

     if (dr_globalMe .eq. MASTER_PE) &
     write(*,*) '###############################################################################'

     ! Compute next step dt:
     ! backup needed old
     dr_dtOld = dr_dt

     ! calculate new
     call Driver_computeDt(dr_nbegin,  dr_nstep,      &
                           dr_simTime, dr_dtOld, dr_dtNew)
     ! store new
     dr_dt = dr_dtNew
     
     call Timers_start("io")
     call IO_output(dr_simTime,dr_dt,dr_nstep+1,dr_nbegin,endRun)
     call Timers_stop("io")

     !!KPD KPD
     !call gr_sbFinalize()

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

  !KPD KPD
  call gr_sbFinalize()

  !KPD KPD
  deallocate(gr_sbBodyInfo)

  call Timers_stop("evolution")
  call Logfile_stamp( 'Exiting evolution loop' , '[Driver_evolveFlash]')
  if(.NOT.endRun) call IO_outputFinal()
  call Timers_getSummary(dr_nstep)
  call Logfile_stamp( "FLASH run complete.", "LOGFILE_END")
  call Logfile_close()

  return
  
end subroutine Driver_evolveFlash



