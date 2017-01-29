!!****if* source/Driver/DriverMain/SolidMechanics_rbc/Driver_evolveFlash
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

  use Driver_data, ONLY: dr_globalMe, dr_globalNumProcs, dr_nbegin,       &
                         dr_nend, dr_dt, dr_wallClockTimeLimit, &
                         dr_tmax, dr_simTime, dr_redshift,      &
                         dr_nstep, dr_dtOld, dr_dtNew,          &
                         dr_restart, dr_elapsedWCTime, dr_globalComm
  use IncompNS_interface, ONLY : IncompNS
  use Driver_interface, ONLY : Driver_sourceTerms, Driver_computeDt, &
       Driver_getElapsedWCTime, Driver_getNstep

  use Logfile_interface,ONLY : Logfile_stamp, Logfile_close
  use Timers_interface, ONLY : Timers_start, Timers_stop, &
                               Timers_getSummary
  use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                                Grid_getListOfBlocks, &
                                Grid_updateRefinement, &
                                Grid_fillGuardCells,  &
                                Grid_getBlkCenterCoords,&
                                Grid_getBlkPtr,   &
                                Grid_releaseBlkPtr, &
                                Grid_getBlkIndexLimits
  use IO_interface,      ONLY : IO_output,IO_outputFinal
  use IO_data , ONLY : io_plotFileNumber, IO_plotFileIntervalStep
  use Grid_data,ONLY : gr_globalDomain,gr_meshNumProcs
  use sm_assemble_interface, ONLY : sm_assemble_IntForce_rbc
  use SolidMechanics_interface, ONLY : SolidMechanics
  use sm_Misc_interface, only: sm_get_NumBodies
  use sm_iointerface, only: sm_ioWriteSnapshot
  use SolidMechanics_data, ONLY :sm_BodyInfo ,sm_meshMe, sm_meshComm,sm_structure
  use SolidMechanics_rbc_data, Only :rbcplotOutputInterval
  use ImBound_interface, ONLY : ImBound
  use gr_sbData , ONLY : gr_sbBodyInfo
  !use Particles_interface, ONLY : Particles_advance, Particles_dump
  !use gr_ptVPData, ONLY : gr_ptVPBndBox, gr_ptVPDeltas
  !use Particles_data, ONLY: particles, pt_numLocal,pt_maxPerProc

  implicit none
   
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
#include "SolidMechanics.h"
#include "Imbound.h"  

  integer :: localNumBlocks  
  integer :: blockCount
  integer :: blockList(MAXBLOCKS)
  integer :: sweepDummy
  integer :: ibd, Numbodies, nstep
  integer :: TA(2),count_rate
  character(len=6) :: index_count,index_proc
  integer :: i,j,k,lb,blockID,ierr
  real :: xi,yi,zi,coords(MDIM)
  character(len=6) :: index_name,index_ibd
  integer count, firstfileflag

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, pointer, dimension(:,:,:,:) :: faceData_low,faceData_high
  integer :: blk_low,blk_high
  real :: aux

  ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(3,2) :: strBuff
  character(len=15) :: numToStr
  real :: ET
  logical :: endRun,restart=.false.  
  type(sm_structure), pointer:: body
  
  endRun = .false.

  
  ! Initial Timestep:
  ! backup needed old
  dr_dtOld = dr_dt
  write(*,*)'dr_dt=',dr_dt
  ! calculate new
  call Driver_computeDt(dr_nbegin,  dr_nstep,      &
       dr_simTime, dr_dtOld, dr_dtNew)
  ! store new
  dr_dt = dr_dtNew
  write(*,*) 'dr_dt=',dr_dt  

  !! Write To Tecplot:
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  count = io_plotFileNumber-1
  firstfileflag = 0  
  !call ins_ioWriteTecplot_grid(dr_globalMe,dr_simtime,dr_dt,dr_nstep,count, &
  !                             0.0,blockList,blockCount,firstfileflag)

  !                  Geomtery written to tecplot file geom.00
  !------------------------------------------------------------
  ! First step (t=0) force calc. 
  call sm_get_NumBodies(NumBodies)
  do ibd=1,numbodies
     call sm_assemble_IntForce_rbc (ibd,SM_IOPT_NSTEP)
     body => sm_bodyinfo(ibd)
     call sm_ioWrite_rbc(ibd,WRITEHDF5)
     
  end do
  
  call Logfile_stamp( 'Entering evolution loop' , '[Driver_evolveFlash]')
  call Timers_start("evolution")
  
  

  loop_integ: do dr_nstep = dr_nBegin, dr_nend
     
     call SYSTEM_CLOCK(TA(1),count_rate)
     
     
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

     dr_simTime = dr_simTime + dr_dt

     call Timers_start("SolidMechanics_advance")
     call SolidMechanics(SM_ADVANCE,restart)
     call ImBound(blockCount, blockList, dr_dt, FORCE_FLOW)
     call ImBound(blockCount, blockList, dr_dt, COMPUTE_FORCES)
     call Timers_stop("SolidMechanics_advance")     


     ! Fill Guardcell for velocities
     call Grid_fillGuardCells(CENTER_FACES,ALLDIR)
     !----
     !- End Physics Sequence
     !--------------------------------------------------------------------
     call Timers_start("Grid_updateRefinement")
     call Grid_updateRefinement( dr_nstep, dr_simTime)
     call Timers_stop("Grid_updateRefinement")
     
     ! HACK: write out a snapshot for the bodies
     call Driver_getNStep(nstep)
     if( mod( nstep, IO_plotFileIntervalStep ) == 0 ) then
        
        ! Write to Snapshot
        call sm_get_NumBodies(NumBodies)
        do ibd = 1,NumBodies
           if (sm_meshMe .eq. sm_BodyInfo(ibd)%BodyMaster)then
              call sm_ioWriteSnapshot(ibd,nstep)
           endif
        end do
         !!! Write to Tecplot - Put in IO file:
        count = count + 1
        write(index_name,"(I6.6)") count
        if (sm_meshMe .eq. MASTER_PE) then
          open(unit=113,file='./IOData/geo.'//index_name//'.plt',form='formatted')
          !write(113,'(A,G12.5,A)')'TEXT X=75,Y=5,F=HELV-BOLD,C=RED,T=" T =',dr_simTime,'"'
          close(113)
        endif
        call MPI_BARRIER(sm_meshComm,ierr)

        do ibd = 1,NumBodies
        if (sm_meshMe .eq. sm_BodyInfo(ibd)%BodyMaster)then
   
           open(unit=113,file='./IOData/geo.'//index_name//'.plt',form='formatted', &
                status='old',position='append')

           write(index_ibd,"(I6.6)") ibd  

           write(113,'(A)') 'VARIABLES = "X", "Y", "Z"'
           write(113,*)                                                &
           'ZONE T=Body'//index_ibd//', N=',                           &
           sm_BodyInfo(ibd)%nnp,', E=',sm_BodyInfo(ibd)%ws_nel,        &
           ', DATAPACKING = POINT, ZONETYPE = FETRIANGLE'
           write(113,*)'DT = (SINGLE SINGLE SINGLE)'
           do i = 1,sm_BodyInfo(ibd)%nnp
              xi= sm_BodyInfo(ibd)%qn(3*i-2)
              yi= sm_BodyInfo(ibd)%qn(3*i-1)
              zi= sm_BodyInfo(ibd)%qn(3*i)   
              write(113,'(3F12.8)') xi,yi,zi
           enddo
           do i=1,sm_BodyInfo(ibd)%ws_nel
              write(113,'(3(I8))')                       &
                    sm_BodyInfo(ibd)%ws_IEN(1,i),        &
                    sm_BodyInfo(ibd)%ws_IEN(2,i),        &
                    sm_BodyInfo(ibd)%ws_IEN(3,i)
           enddo

           close(113)

        endif
        call MPI_BARRIER(sm_meshComm,ierr)
        enddo
        !!!! -----------


        !! Write Grid To Tecplot:
        call Grid_getListOfBlocks(LEAF,blockList,blockCount)
        !call ins_ioWriteTecplot_grid(dr_globalMe,dr_simtime,dr_dt,dr_nstep,count, &
         !                            0.0,blockList,blockCount,firstfileflag)

     end if

     
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
     
     if((dr_nstep==1).or.(MOD(dr_nstep,rbcplotOutputInterval)==0))  then
        do ibd=1,Numbodies
           write(*,*) 'Writing outputfiles at iteration....',dr_nstep
           call Timers_start("WRITEHDF5")
           call sm_ioWrite_rbc(ibd,WRITEHDF5)
           call Timers_stop("WRITEHDF5")
        end do
     endif
     

     
     
     
     if (dr_globalMe == MASTER_PE) then
     
        
        !  Write the step time to output
        CALL SYSTEM_CLOCK(TA(2),count_rate)
        ET=REAL(TA(2)-TA(1),8)/count_rate  
        write(*,*)'  '
        write(*,*)'!---------------------------------------------'
        write(*,'(A,I7,A,F12.8,A,F12.8)')'Step=',int(dr_nstep),'  t=',dr_simTime,'  dt=',dr_dt
        write(*,*)'!---------------- '
        write(*,'(A,F12.8,A)')'Total Step Time =',ET, ' sec. '
        write(*,*)'!---------------------------------------------'
        write(*,*)'  '
        if (dr_simTime >= dr_tmax) then
           print *, "exiting: reached max SimTime"
           exit
        endif
        if(endRun) exit
     end if
     
  enddo loop_integ
  
  call Timers_stop("evolution")
  call Logfile_stamp( 'Exiting evolution loop' , '[Driver_evolveFlash]')
  if(.NOT.endRun) call IO_outputFinal()
  call Timers_getSummary( dr_nstep)
  call Logfile_stamp( "FLASH run complete.", "LOGFILE_END")
  call Logfile_close()
  
  return
  
end subroutine Driver_evolveFlash



!!$  real,save :: Domainsize(3)
!!$  ! Calculate the domain size to be used later in case of PERIODIC BC.
!!$  domainsize(IAXIS)=gr_globalDomain(HIGH,IAXIS)- gr_globalDomain(LOW,IAXIS)
!!$  domainsize(JAXIS)=gr_globalDomain(HIGH,JAXIS)- gr_globalDomain(LOW,JAXIS)
!!$#if NDIM==3
!!$  domainsize(KAXIS)=gr_globalDomain(HIGH,KAXIS)- gr_globalDomain(LOW,KAXIS)
!!$#endif
!!$  if (dr_globalMe==MASTER_PE) write(*,*)'Domain dimensions=',domainsize(IAXIS:KAXIS)
!!$  
!!$  !--------------------
!!$  ! write the geomtry to tecplot
!!$  !----------------------
!!$  call Grid_getLocalNumBlks(localNumBlocks)
!!$  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
!!$  call int2char(dr_globalMe,index_proc)
!!$  open(unit=114,file='./IOData/geom'//'.'//index_proc//'.plt',form='formatted')
!!$  write(gridname,'("geom",I5.5)') int(dr_nstep)
!!$  write(114,*)'TITLE = "Domain grid "'
!!$  write(114,*)'VARIABLES = "X", "Y", "Z"'  
!!$  write(114,'("ZONE T=",A," DATAPACKING=POINT, N=",I0,", E=",I0,", ZONETYPE=FEBRICK")') TRIM(gridname),8*blockCount, blockCount 
!!$  do i=1,blockCount
!!$     blkID=blocklist(i);
!!$     call Grid_getBlkBoundBox(blkID,gr_ptVPBndBox)
!!$     write(114,*) gr_ptVPBndBox(LOW,IAXIS),gr_ptVPBndBox(LOW,JAXIS),gr_ptVPBndBox(LOW,KAXIS)
!!$     write(114,*) gr_ptVPBndBox(HIGH,IAXIS),gr_ptVPBndBox(LOW,JAXIS),gr_ptVPBndBox(LOW,KAXIS)
!!$     write(114,*) gr_ptVPBndBox(HIGH,IAXIS),gr_ptVPBndBox(HIGH,JAXIS),gr_ptVPBndBox(LOW,KAXIS)
!!$     write(114,*) gr_ptVPBndBox(LOW,IAXIS),gr_ptVPBndBox(HIGH,JAXIS),gr_ptVPBndBox(LOW,KAXIS)
!!$     write(114,*) gr_ptVPBndBox(LOW,IAXIS),gr_ptVPBndBox(LOW,JAXIS),gr_ptVPBndBox(HIGH,KAXIS)
!!$     write(114,*) gr_ptVPBndBox(HIGH,IAXIS),gr_ptVPBndBox(LOW,JAXIS),gr_ptVPBndBox(HIGH,KAXIS)
!!$     write(114,*) gr_ptVPBndBox(HIGH,IAXIS),gr_ptVPBndBox(HIGH,JAXIS),gr_ptVPBndBox(HIGH,KAXIS)
!!$     write(114,*) gr_ptVPBndBox(LOW,IAXIS),gr_ptVPBndBox(HIGH,JAXIS),gr_ptVPBndBox(HIGH,KAXIS)
!!$  end do
!!$  write(114,*) '  '
!!$  do i=1,blockCount
!!$     write(114,'(8(I5))') (/(j,j=1+(i-1)*8,i*8)/)
!!$  end do
!!$  close(114)


