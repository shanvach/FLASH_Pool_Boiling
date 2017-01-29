!!****if* source/Simulation/SimulationMain/VParticles_DPD/Driver_evolveFlash
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
!!  This is a modification of the standard Driver_evolveFlash for the unitTest for 
!!     the Particles unit.  Instead of using Hydro sweeps in xyz/zyx, the
!!     velocities are generated with a random perturbation overlaid on a 
!!     constant flow in the x-direction.
!!
!! NOTES
!!
!!
!!
!!***

! Debugging flags
#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif
!#define DEBUG_VPARTICLES



subroutine Driver_evolveFlash()

  use Driver_data, ONLY: dr_meshNumProcs,dr_globalMe, dr_nbegin, &
       dr_nend, dr_dt, dr_wallClockTimeLimit, &
       dr_tmax, dr_simTime, dr_redshift, &
       dr_nstep, dr_dtOld, dr_dtNew, dr_nbegin, dr_restart
  use Driver_interface, ONLY : Driver_computeDt,Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stamp, Logfile_close
  use Timers_interface, ONLY : Timers_start, Timers_stop, &
    Timers_getSummary
  use Particles_interface, ONLY : Particles_advance, &
    Particles_unitTest
  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
    Grid_getListOfBlocks, Grid_updateRefinement
  use gr_ptVPData, ONLY : gr_ptVPBndBox, gr_ptVPDeltas
  use IO_interface, ONLY : IO_output, IO_outputFinal
  use Particles_data, ONLY: particles, pt_numLocal,pt_maxPerProc
  use Simulation_data, ONLY: domainsize,sqrt_dt,pt_NumPart
  use Grid_data,ONLY : gr_globalDomain,gr_meshNumProcs
  implicit none
 
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  !------- Local variables -------------------
  
  integer, parameter :: stepsPerAdvance = 2
  integer, parameter :: fileUnit = 2
  integer,dimension(3) :: opind,npind,ovind,nvind,intvind,ofind,nfind
  integer,dimension(4) :: prNum
  integer :: blockCount
  integer :: blockList(MAXBLOCKS)
  integer :: localNumBlocks, i,j
  integer :: temp
  integer :: dummy5(5) !dummy val to pass into Particles_computeDt, this arg not needed here
  integer :: kk, blockID
  integer :: count,blkID,p_count
  integer :: TA(2),count_rate
  integer ::  comm, ierr, total_pt_local
  integer, save :: firstfileflag
  logical ::  perfect = .true.
  logical :: endRun
  logical :: test1(3)

  character(len=6) :: index_count,index_proc
  character(len=100) :: fileName
  character(len=20) :: gridname
  character(len=MAX_STRING_LENGTH), dimension(3,2) :: strBuff
  character(len=15) :: numToStr

  real,dimension(NDIM,pt_maxPerProc):: pos,npos,v,vnew,f,fnew,v_telda 
  real :: dt
  real :: ET
!------ End of local variables -------------
  
  
  ! stays true if no errors are found
  perfect = .true.
  call Logfile_stamp( 'Entering evolution loop' , '[DRIVER_GLOBAL]')
  call Timers_start("evolution")
  !----
  fileName='DPD_particles_'
  count=0
  firstfileflag = 1 ! io_plotFileNumber

  ! Calculate the domain size to be used later in case of PERIODIC BC.
  domainsize(IAXIS)=gr_globalDomain(HIGH,IAXIS)- gr_globalDomain(LOW,IAXIS)
  domainsize(JAXIS)=gr_globalDomain(HIGH,JAXIS)- gr_globalDomain(LOW,JAXIS)
#if NDIM==3
  domainsize(KAXIS)=gr_globalDomain(HIGH,KAXIS)- gr_globalDomain(LOW,KAXIS)
#endif
  if (dr_globalMe==MASTER_PE) write(*,*)'Domain dimensions=',domainsize(IAXIS:KAXIS)
  
  !--------------------
  ! write the geomtry to tecplot
  !----------------------
  call Grid_getLocalNumBlks(localNumBlocks)
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  call int2char(dr_globalMe,index_proc)

  p_count=pt_numLocal;
  call IO_OutputDPDtoParaview (filename,dr_globalMe,dr_simtime,dr_dt,dr_nstep,count, &
     blockList,blockCount,firstfileflag,particles,p_count )

  firstfileflag = 0
  open(unit=114,file='./IOData/geom'//'.'//index_proc//'.plt',form='formatted')
  write(gridname,'("geom",I5.5)') int(dr_nstep)
  write(114,*)'TITLE = "Domain grid "'
  write(114,*)'VARIABLES = "X", "Y", "Z"'  
  write(114,'("ZONE T=",A," DATAPACKING=POINT, N=",I0,", E=",I0,", ZONETYPE=FEBRICK")') TRIM(gridname),8*blockCount, blockCount 

  do i=1,blockCount
     blkID=blocklist(i);
     call Grid_getBlkBoundBox(blkID,gr_ptVPBndBox)
     write(114,*) gr_ptVPBndBox(LOW,IAXIS),gr_ptVPBndBox(LOW,JAXIS),gr_ptVPBndBox(LOW,KAXIS)
     write(114,*) gr_ptVPBndBox(HIGH,IAXIS),gr_ptVPBndBox(LOW,JAXIS),gr_ptVPBndBox(LOW,KAXIS)
     write(114,*) gr_ptVPBndBox(HIGH,IAXIS),gr_ptVPBndBox(HIGH,JAXIS),gr_ptVPBndBox(LOW,KAXIS)
     write(114,*) gr_ptVPBndBox(LOW,IAXIS),gr_ptVPBndBox(HIGH,JAXIS),gr_ptVPBndBox(LOW,KAXIS)
     write(114,*) gr_ptVPBndBox(LOW,IAXIS),gr_ptVPBndBox(LOW,JAXIS),gr_ptVPBndBox(HIGH,KAXIS)
     write(114,*) gr_ptVPBndBox(HIGH,IAXIS),gr_ptVPBndBox(LOW,JAXIS),gr_ptVPBndBox(HIGH,KAXIS)
     write(114,*) gr_ptVPBndBox(HIGH,IAXIS),gr_ptVPBndBox(HIGH,JAXIS),gr_ptVPBndBox(HIGH,KAXIS)
     write(114,*) gr_ptVPBndBox(LOW,IAXIS),gr_ptVPBndBox(HIGH,JAXIS),gr_ptVPBndBox(HIGH,KAXIS)
  end do
  write(114,*) '  '
  do i=1,blockCount
     write(114,'(8(I5))') (/(j,j=1+(i-1)*8,i*8)/)
  end do
  close(114)
  
  !                  Geomtery written to tecplot file geom.00
  !------------------------------------------------------------
  
  !-------------------------------------------------------------
  !                     First step calculations 

  dt=dr_dt;
  sqrt_dt=sqrt(dr_dt);

  !This is a call to particle advace to make sure that all the particles are copied in there correct processors
!!$  write(*,*) 'pt_numLocal',pt_numLocal,'p_count',p_count
!!$    do i=1,pt_numLocal
!!$     write(*,*) particles(POSX_PART_PROP:POSZ_PART_PROP,i)
!!$  end do

  write(*,*) 'BEFORE PARTICLES ADVANCE'

  call Particles_advance(dr_dt,dr_dt)
  p_count=pt_numLocal;
  write(*,*) 'AFTER PARTICLES ADVANCE'

!!$  do i=1,p_count
!!$     write(*,*) particles(POSX_PART_PROP:POSZ_PART_PROP,i)
!!$  end do

!!$  ! Update the simTime
  dr_simTime=dr_simTime+dr_dt;
  
  !----------------------------------------------------------------
  !                      The iterations loop
  
  loop_integ: do dr_nstep = dr_nbegin, dr_nend
     count=count+1
     call SYSTEM_CLOCK(TA(1),count_rate)

     if (dr_globalMe == MASTER_PE) then
        write(*,*)'  '
        write(*,*)'!---------------- '
        write(*,'(A,I7,A,F12.8,A,F12.8)')'Step=',int(dr_nstep),'  t=',dr_simTime,'  dt=',dt
        write(*,*)'!---------------- '
        write(*,*)' '
        call Grid_getLocalNumBlks(localNumBlocks)
        call Grid_getListOfBlocks(LEAF,blockList,blockCount)
        
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
          
     p_count=pt_numLocal;
     
     write(*,*)'Local number of particles on processor ',dr_globalMe, p_count

     call Particles_advance(dr_dt,dr_dt) 
     !print *,'done with Particles Advance again'
     
     if ( dr_meshNumProcs>1) then
        total_pt_local=0;
        call MPI_REDUCE(pt_numlocal, total_pt_local,CONSTANT_ONE, & 
             MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        
        if (dr_globalMe==MASTER_PE) then 
           write(*,*)'Total no. of particles= ',total_pt_local
           if (total_pt_local/=pt_NumPart) stop
        end if
     end if

     
     !write(*,*)particles(:,1)
     call IO_writeDpd(particles,pt_numLocal)

     
     !call Timers_start("Grid_updateRefinement")
     !call Grid_updateRefinement( dr_nstep, dr_simTime)
     !call Timers_stop("Grid_updateRefinement")
     
     call IO_output( dr_simTime, dr_dt, dr_nstep+1, dr_nbegin, endRun)
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
        end if
        exit
     end if
     
     if (dr_simTime >  dr_wallClockTimeLimit) then
        if(dr_globalMe == MASTER_PE) then
           print *, "exiting: reached max wall clock time"
        end if
        exit
     end if
     
     
     dr_dtOld = dr_dt                                  ! backup needed old 
     
     ! calculate new timestep, needed for Particles_advance
     !! don't use Particles_computeDt because this only works with the
     !! local number of blocks.  Driver_computeDt takes care of
     !! MPI calls. need to loop over all blocks
     
     
     call Driver_computeDt( dr_nBegin, dr_nStep, &
          dr_simTime, dr_dtOld, dr_dtNew)
     !   write(*,*)'New sim time=',dr_simTime
     
     dr_dt = dr_dtNew                                    ! store new
     dr_simTime=dr_simTime+dr_dt;
     sqrt_dt=sqrt(dr_dt);

     !  Write the step time to output
     CALL SYSTEM_CLOCK(TA(2),count_rate)
     ET=REAL(TA(2)-TA(1),8)/count_rate
     if (dr_globalMe.eq. MASTER_PE)  then 
        write(*,*)'  '
        write(*,*)'!---------------------------------------------'
        write(*,'(A,F12.8)')'Total Step Time =',ET
        write(*,*)'!---------------------------------------------'
        write(*,*)'  '
     end if
  end do loop_integ
  
  call Timers_stop("evolution")
  
  call Logfile_stamp( 'Exiting evolution loop' , '[DRIVER_GLOBAL]')
  
  if(.NOT.endRun) call IO_outputFinal( )
     
     call Timers_getSummary(dr_nstep)
     
    
     !finish unit test write out file
     
     if (perfect) then
        write(fileUnit,'("Particles unitTest PASSED!  all results conformed with expected values.")')
        !write(*,'("Particles unitTest PASSED!  all results conformed with expected values.")')
        call Logfile_stamp( "Particles unitTest PASSED!")
     else
        write(fileUnit,'("Particles unitTest FAILED!")')
        write(*,'("Particles unitTest FAILED!")')
        call Logfile_stamp( "Particles unitTest FAILED!")
     endif
     
     close(fileUnit)
     
     call Logfile_stamp( "FLASH run complete.", "LOGFILE_END")
     call Logfile_close()
     
     if (dr_globalMe==MASTER_PE) write(*,*)"FLASH run complete."
     return
     
   end subroutine Driver_evolveFlash
   
   
   
   Subroutine int2char(i,strng)
     
     integer i
     character (6) strng
     
     integer k, val, valaux
     real*8 val2
     
     valaux=0 
    strng = '000000'
    
    
    do k = 6,1,-1
       
       val2 = (i-valaux) / (10**(k-1))
       val = floor(val2) 
       
       valaux = valaux + val*(10**(k-1))
       
       !        write(*,*) 7-k,val,valaux
       
         if (val .GE. 1) then
            
            select case (val)
               
            case (1)
               
               strng(7-k:7-k) = "1"
               
            case (2)
               
               strng(7-k:7-k) = '2'
               
               
            case (3)
               
               strng(7-k:7-k) = '3'
               
            case (4)
               
               strng(7-k:7-k) = '4'
               
            case (5)
               
               strng(7-k:7-k) = '5'
               
            case (6)
               
               strng(7-k:7-k) = '6'
               
            case (7)
               
               strng(7-k:7-k) = '7'
               
            case (8)
               
               strng(7-k:7-k) = '8'
               
            case (9)
               
               strng(7-k:7-k) = '9'
               
            end select
            
         endif

      enddo

      End subroutine int2char


subroutine  MACHINEEPSILON
  use Simulation_data, ONLY: eps
  implicit none
  real :: MACHEPS=1.D0

  MACHEPS = MACHEPS / 2.D0
  !MACHEPS = 1.D0
  do while (1.D0 + MACHEPS / 2.D0 .NE. 1.D0) 
     MACHEPS = MACHEPS / 2.D0
     !  IF ( 1.D0 + MACHEPS / 2.D0 .EQ. 1.D0 ) GOTO 110
  end do
  !GO TO 100
  !110  CONTINUE
  eps= MACHEPS;
  print*, MACHEPS
  return
end subroutine MACHINEEPSILON
