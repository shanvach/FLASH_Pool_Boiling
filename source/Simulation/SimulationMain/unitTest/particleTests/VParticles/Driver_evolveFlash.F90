!!****if* source/Simulation/SimulationMain/unitTest/particleTests/VParticles/Driver_evolveFlash
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


#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif

subroutine Driver_evolveFlash()

  use Driver_data, ONLY: dr_globalMe, dr_nbegin, &
       dr_nend, dr_dt, dr_wallClockTimeLimit, &
       dr_tmax, dr_simTime, dr_redshift, &
       dr_nstep, dr_dtOld, dr_dtNew, dr_nbegin, dr_restart
  use Driver_interface, ONLY : Driver_computeDt
  use Logfile_interface, ONLY : Logfile_stamp, Logfile_close
  use Timers_interface, ONLY : Timers_start, Timers_stop, &
    Timers_getSummary
  use Particles_interface, ONLY : Particles_advance, &
    Particles_unitTest
  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
    Grid_getListOfBlocks, Grid_updateRefinement
  use gr_ptVPData, ONLY : gr_ptVPBndBox, gr_ptVPDeltas
  use IO_interface, ONLY : IO_output, IO_outputFinal
  use Particles_data, ONLY: particles, pt_numLocal
  implicit none

#include "constants.h"
#include "Flash.h"


  integer   :: localNumBlocks, i

  integer, parameter :: stepsPerAdvance = 2

  integer :: blockCount
  integer :: blockList(MAXBLOCKS)

  
  ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(3,2) :: strBuff
  character(len=15) :: numToStr


  logical ::  perfect = .true.

  character(len=27) :: fileName
  integer, parameter        :: fileUnit = 2
  integer,dimension(4) :: prNum
  integer :: temp
  integer :: dummy5(5) !dummy val to pass into Particles_computeDt, this arg not needed here
  integer :: kk, blockID
  logical :: endRun
  character(len=6) :: index_count,index_proc

  integer :: count,blkID

  ! stays true if no errors are found
  perfect = .true.
  
  temp = dr_globalMe
  !write(*,*)'dr_globalMe= ',dr_globalMe
  
  do i = 1,4
     prNum(i)= mod(temp,10)
     temp = temp/10
  end do


  filename = "unitTest_"//char(48+prNum(4))//char(48+prNum(3))//&
                                 char(48+prNum(2))//char(48+prNum(1))
  !write(*,*)'filename= ',filename  
  open(fileUnit,file=fileName)
  write(fileUnit,'("P",I0)') dr_globalMe


  call Logfile_stamp( 'Entering evolution loop' , '[DRIVER_GLOBAL]')

  call Timers_start("evolution")
  count=0

#ifdef DEBUG_VPARTICLES
     call int2char(count,index_count)
     call int2char(dr_globalMe,index_proc)
     open(unit=113,file='./IOData/part'//index_count//'.'//index_proc//'.plt',form='formatted')
     !write(113,'(A,G12.5,A)')'TEXT X=75,Y=5,F=HELV-BOLD,C=RED,T=" T =',dr_simTime,'"'

       
        write(113,'(A)') 'VARIABLES = "X" , "Y"'
        write(113,'(A,I8,A)')'ZONE I=',pt_numLocal,',DATAPACKING = POINT' 

        do i = 1,pt_numLocal
           write(113,'(2E16.8)') particles(POSX_PART_PROP,i),particles(POSY_PART_PROP,i)
        enddo

     close(113)


     call Grid_getLocalNumBlks(localNumBlocks)
     call Grid_getListOfBlocks(LEAF,blockList,blockCount)
     print *,' block count is ',blockCount
     print *,' LeaF blocks are ',blocklist(1:blockCount)
     do i=1,blockCount
       blkID=blocklist(i);
       call Grid_getDeltas(blkID,gr_ptVPDeltas)
       call Grid_getBlkBoundBox(blkID,gr_ptVPBndBox)
       print *,' Proc ',dr_globalMe,' Block ',blocklist(i)
       print *,'Xmin=',gr_ptVPBndBox(LOW,IAXIS),'Xmax=',gr_ptVPBndBox(HIGH,IAXIS)
       print *,'Ymin=',gr_ptVPBndBox(LOW,JAXIS),'Ymax=',gr_ptVPBndBox(HIGH,JAXIS)
     enddo

#endif   

     do dr_nstep = dr_nbegin, dr_nend
        !!     print *,' at step ',dr_nstep
        
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
     
     call Particles_advance(dr_dt,dr_dt) 
     !print *,'done with Particles Advance again'
     
     call Particles_clean(particles,pt_numLocal)
     !print *,'Done with particles cleaning'
     !print *,'This is the new local particles number ',pt_numLocal

#ifdef DENUG_VPARTICLES
     count=count+1
     call int2char(count,index_count)
     call int2char(dr_globalMe,index_proc)
     open(unit=113,file='./IOData/part'//index_count//'.'//index_proc//'.plt',form='formatted')
     !write(113,'(A,G12.5,A)')'TEXT X=75,Y=5,F=HELV-BOLD,C=RED,T=" T =',dr_simTime,'"'

       
        write(113,'(A)') 'VARIABLES = "X" , "Y"'
        write(113,'(A,I8,A)')'ZONE I=',pt_numLocal,',DATAPACKING = POINT' 

        do i = 1,pt_numLocal
           write(113,'(2E16.8)') particles(POSX_PART_PROP,i),particles(POSY_PART_PROP,i)
        enddo

     close(113)

#endif
     !! Test whether particles are living outside their local blocks
     call Particles_unitTest(fileUnit,perfect)
!!     print *,'done with Particles_unitTest'

     call Timers_start("Grid_updateRefinement")
     call Grid_updateRefinement( dr_nstep, dr_simTime)
     call Timers_stop("Grid_updateRefinement")

     call IO_output( dr_simTime, &
          dr_dt, dr_nstep+1, dr_nbegin, endRun)
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
     

     dr_dtOld = dr_dt                     ! backup needed old 

     ! calculate new timestep, needed for Particles_advance
     !! don't use Particles_computeDt because this only works with the
     !! local number of blocks.  Driver_computeDt takes care of
     !! MPI calls
     ! need to loop over all blocks


!!     print *, "before Driver_computeDT"
     call Driver_computeDt( dr_nBegin, dr_nStep, &
                    dr_simTime, dr_dtOld, dr_dtNew)

!!     print *,'back from  compute Particles_computeDt'

     dr_dt = dr_dtNew                                    ! store new



  enddo

  call Timers_stop("evolution")

  call Logfile_stamp( 'Exiting evolution loop' , '[DRIVER_GLOBAL]')

  if(.NOT.endRun) call IO_outputFinal( )

  call Timers_getSummary( dr_nstep)




  !finish unit test write out file

  if (perfect) then
    write(fileUnit,'("Particles unitTest PASSED!  all results conformed with expected values.")')
    write(*,'("Particles unitTest PASSED!  all results conformed with expected values.")')
    call Logfile_stamp( "Particles unitTest PASSED!")
  else
    write(fileUnit,'("Particles unitTest FAILED!")')
    write(*,'("Particles unitTest FAILED!")')
    call Logfile_stamp( "Particles unitTest FAILED!")
  endif

  close(fileUnit)

  call Logfile_stamp( "FLASH run complete.", "LOGFILE_END")
  call Logfile_close()


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

