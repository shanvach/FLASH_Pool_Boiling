!!****if* source/Simulation/SimulationMain/SolidMechanics/First_test_Rigid/Driver_evolveFlash
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
                         dr_restart, dr_elapsedWCTime,dr_globalComm

  use IncompNS_interface, ONLY : IncompNS
  use Driver_interface, ONLY : Driver_sourceTerms, Driver_computeDt, &
       Driver_getElapsedWCTime, Driver_getNStep
  use Logfile_interface,ONLY : Logfile_stamp, Logfile_close
  use Timers_interface, ONLY : Timers_start, Timers_stop, &
                               Timers_getSummary
  use Particles_interface, ONLY : Particles_advance, Particles_dump
  use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                                Grid_getListOfBlocks, &
                                Grid_updateRefinement,&
                                Grid_fillGuardCells,  &
                                Grid_getBlkCenterCoords,&
                                Grid_getBlkPtr,   &
                                Grid_releaseBlkPtr, &
                                Grid_getBlkIndexLimits
  use IO_interface,      ONLY : IO_output,IO_outputFinal
  use SolidMechanics_interface, ONLY: SolidMechanics
  use sm_Misc_interface, only: sm_get_NumBodies
  use sm_iointerface, only: sm_ioWriteSnapshot


  use SolidMechanics_Data, only : sm_meshMe,sm_meshComm,sm_BodyInfo
  use IO_data , ONLY : io_plotFileNumber, IO_plotFileIntervalStep

  use ImBound_interface, ONLY : ImBound

  use gr_sbData , ONLY : gr_sbBodyInfo

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
#include "SolidMechanics.h"
#include "ImBound.h"

  integer   :: localNumBlocks

  integer :: blockCount
  integer :: blockList(MAXBLOCKS)
  integer :: sweepDummy

  integer :: ibd, NumBodies, nstep
  
  ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(3,2) :: strBuff
  character(len=15) :: numToStr

  logical :: endRun  

  integer :: i,j,k,lb,blockID,ierr
  real :: xi,yi,zi,coords(MDIM)
  character(len=6) :: index_name,index_ibd

  integer count, firstfileflag

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, pointer, dimension(:,:,:,:) :: faceData_low,faceData_high
  integer :: blk_low,blk_high
  real :: aux

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

  !! Write To Tecplot:
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  count = io_plotFileNumber-1
  firstfileflag = 0  
  call ins_ioWriteTecplot_grid(dr_globalMe,dr_simtime,dr_dt,dr_nstep,count, &
                               0.0,blockList,blockCount,firstfileflag)
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

!!$     ! Coordinates of Blocks:
!!$     do lb=1,blockCount
!!$        blockID = blockList(lb)
!!$        call Grid_getBlkCenterCoords(blockID,coords)
!!$        !print*, blockCount, blockID, coords(1:MDIM)
!!$     enddo


!!$     ! Spit differences on block boundary velocities:
!!$     ! Get Blocks internal limits indexes:
!!$     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC) 
!!$     ! U 
!!$     blk_high = 1
!!$     blk_low  = 2
!!$     call Grid_getBlkptr(blk_low,faceData_low  ,FACEX)
!!$     call Grid_getBlkptr(blk_high,faceData_high,FACEX)
!!$     write(*,*) 'Diff 1h 2L='
!!$     do k=1,blkLimitsGC(HIGH,KAXIS) 
!!$        do j=1,blkLimitsGC(HIGH,KAXIS)
!!$           aux = faceData_high(VELC_FACE_VAR,blkLimits(HIGH,IAXIS)+1,j,k)- &
!!$                 faceData_low(VELC_FACE_VAR,blkLimits(LOW,IAXIS),j,k)
!!$           if (abs(aux) .gt. 1.e-10) write(*,*) j,k,aux
!!$        enddo
!!$     enddo
!!$     call Grid_releaseBlkptr(blk_low,faceData_low  ,FACEX)
!!$     call Grid_releaseBlkptr(blk_high,faceData_high,FACEX)



     !--------------------------------------------------------------------
     !- Start Physics Sequence
     !----
#ifdef DEBUG_DRIVER
     print*, 'going into SolidMechanics'
#endif
     dr_simTime = dr_simTime + dr_dt

     call Timers_start("SolidMechanics")
     call SolidMechanics(SM_ADVANCE1DT)

     call ImBound(blockCount, blockList, dr_dt, FORCE_FLOW)
     call ImBound(blockCount, blockList, dr_dt, COMPUTE_FORCES)
     call Timers_stop("SolidMechanics")

     ! Fill Guardcell for velocities
     call Grid_fillGuardCells(CENTER_FACES,ALLDIR)

     ! Check what we have at the end of forcing in velocities:
     ! Spit differences on block boundary velocities:
!!$     ! U 
!!$     blk_high = 1
!!$     blk_low  = 2
!!$     call Grid_getBlkptr(blk_low,faceData_low  ,FACEX)
!!$     call Grid_getBlkptr(blk_high,faceData_high,FACEX)
!!$     write(*,*) 'Diff 1h 2L After forcing='
!!$     do k=1,blkLimitsGC(HIGH,KAXIS) 
!!$        do j=1,blkLimitsGC(HIGH,KAXIS)
!!$           aux = faceData_high(VELC_FACE_VAR,blkLimits(HIGH,IAXIS)+1,j,k)- &
!!$                 faceData_low(VELC_FACE_VAR,blkLimits(LOW,IAXIS),j,k)
!!$           if (abs(aux) .gt. 1.e-10) write(*,*) j,k,aux
!!$        enddo
!!$     enddo
!!$     call Grid_releaseBlkptr(blk_low,faceData_low  ,FACEX)
!!$     call Grid_releaseBlkptr(blk_high,faceData_high,FACEX)

!!$     call sm_get_NumBodies(NumBodies)
!!$     do ibd = 1,NumBodies
!!$        do i=1,gr_sbBodyInfo(ibd)%totalpart
!!$           write(*,*) i,gr_sbBodyInfo(ibd)%particles(PRES_PART_PROP,i)
!!$        enddo
!!$     enddo

#ifdef DEBUG_DRIVER
  print*, 'return from SolidMechanics timestep'
#endif

!!$     call Timers_start("Particles_advance")
!!$     call Particles_advance(dr_dtOld, dr_dt)
!!$#ifdef DEBUG_DRIVER
!!$     print*, 'return from Particles_advance '
!!$#endif
!!$     call Timers_stop("Particles_advance")     

     !----
     !- End Physics Sequence
     !--------------------------------------------------------------------

     call Timers_start("Grid_updateRefinement")
     call Grid_updateRefinement( dr_nstep, dr_simTime)
     call Timers_stop("Grid_updateRefinement")
     
     call Timers_start("io")
     call IO_output(dr_simTime,dr_dt,dr_nstep+1,dr_nbegin,endRun)
     call Timers_stop("io")

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
              xi= sm_BodyInfo(ibd)%x(i) + sm_BodyInfo(ibd)%qn(sm_BodyInfo(ibd)%ID(IAXIS,i))
              yi= sm_BodyInfo(ibd)%y(i) + sm_BodyInfo(ibd)%qn(sm_BodyInfo(ibd)%ID(JAXIS,i))
              zi= sm_BodyInfo(ibd)%z(i) + sm_BodyInfo(ibd)%qn(sm_BodyInfo(ibd)%ID(KAXIS,i))
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
        call ins_ioWriteTecplot_grid(dr_globalMe,dr_simtime,dr_dt,dr_nstep,count, &
                                     0.0,blockList,blockCount,firstfileflag)





     end if

     ! backup needed old
     dr_dtOld = dr_dt

     ! calculate new
     call Driver_computeDt( dr_nbegin,  dr_nstep,      &
                            dr_simTime, dr_dtOld, dr_dtNew)
     ! store new
     dr_dt = dr_dtNew

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
  call Logfile_stamp( 'Exiting evolution loop' , '[Driver_evolveFlash]')
  if(.NOT.endRun) call IO_outputFinal()
  call Timers_getSummary(dr_nstep)
  call Logfile_stamp( "FLASH run complete.", "LOGFILE_END")
  call Logfile_close()

  return
  
end subroutine Driver_evolveFlash



