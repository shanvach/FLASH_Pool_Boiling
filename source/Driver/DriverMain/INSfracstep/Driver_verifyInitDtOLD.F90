!!****if* source/Driver/DriverMain/INSfracstep/Driver_verifyInitDt
!!
!! NAME
!!  Driver_verifyInitDt
!!
!! SYNOPSIS
!!  Driver_verifyInitDt()
!!
!! DESCRIPTION
!!
!!  Verifies that the initial dt passes CFL criteria.
!!
!! ARGUMENTS
!!
!!
!! NOTES 
!!
!!  variables that begin with "dr_" like, dr_dt or dr_globalMe
!!  are stored in the data fortran module for the Driver unit, Driver_data.
!!  The "dr_" is meant to indicate that the variable belongs to the Driver Unit.
!!  all other normally named variables i, j, etc are local variables.
!!
!!***

subroutine Driver_verifyInitDt()

  use Driver_data, ONLY : dr_dt, dr_restart, dr_dtOld, dr_dtInit, dr_globalMe

  use Grid_interface, ONLY : Grid_getListOfBlocks, &
                             Grid_getBlkIndexLimits, &
                             Grid_getCellCoords, &
                             Grid_getDeltas, &
                             Grid_getBlkPtr, &
                             Grid_releaseBlkPtr

  use IncompNS_interface, ONLY : IncompNS_computeDt

  implicit none       

#include "Flash.h"
#include "constants.h" 
#include "IncompNS.h"

  include "Flash_mpi.h"
  
  real :: dtCheck  ,dtCFL
  integer :: localNumBlocks

  integer    :: dtMinLoc(5)
  integer :: i, blockID, ierr
  integer, dimension(MAXBLOCKS) :: blockList

  integer :: coordSize
  logical :: gcell = .true.
  real, dimension(MDIM) :: del



  !arrays which hold the starting and ending indicies of a block
  integer,dimension(2,MDIM)::blkLimits,blkLimitsGC

  !!coordinate infomration to be passed into physics  
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData
  integer :: isize,jsize,ksize,istat


  if (.not. dr_restart) then
     ! compute the CFL timestep for the simulation and compare it to the
     ! user specified initial timestep.  Scream loudly if there is a problem.

     !initialize values 
     dtCheck = huge(dtCheck)
     
     call Grid_getListOfBlocks(LEAF,blockList,localNumBlocks)

     do i = 1, localNumBlocks

        !There is some overhead in calling Hydro_computeDt.  Although it is a
        !pain to get the coordinates and solution data before calling the 
        !routine, this is just initialization.  Getting the coordinates inside
        !Hydro_computeDt would be much more costly during the run
        
        
        !!Get the coordinate information for all the
        call Grid_getBlkIndexLimits(blockList(i),blkLimits,blkLimitsGC)
        isize = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
        jsize = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
        ksize = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1


        blockID=blockList(i)

        
        call Grid_getDeltas(blockID, del)

        call Grid_getBlkPtr(blockID, facexData,FACEX)
        call Grid_getBlkPtr(blockID, faceyData,FACEY)
        call Grid_getBlkPtr(blockID, facezData,FACEZ)


        call IncompNS_computeDt(blockID,             & 
                           isize, jsize, ksize,                &
                           del(DIR_X), del(DIR_Y), del(DIR_Z), &
                           blkLimits,blkLimitsGC,              &
                           facexData,faceyData,                &
                           facezData,                          &
                           dtCheck )

        call Grid_releaseBlkPtr(blockID, facexData, FACEX)
        call Grid_releaseBlkPtr(blockID, faceyData, FACEY)
        call Grid_releaseBlkPtr(blockID, facezData, FACEZ) 
   
     end do


     ! find the minimum across all processors, store it in dtCFL on MasterPE
     call MPI_AllReduce(dtCheck, dtCFL, 1, FLASH_REAL, MPI_MIN, &
          MPI_COMM_WORLD, ierr)

     if (dr_dtInit > 0.1*dtCFL) then
        
        if (dr_globalMe .EQ. MASTER_PE) then
           print *, '***********************************************************'
           print *, ' Warning: The initial timestep is too large.'
           print *, '   initial timestep = ', dr_dtInit
           print *, '   CFL timestep     = ', dtCFL
           print *, ' Resetting dt_init to 0.1*dt_cfl.'
           print *, '***********************************************************'
           print *, ' '
        endif
        
        dr_dt = 0.1*dtCFL
        
     else
        
        dr_dt = dr_dtInit
        
     endif
     
     dr_dtOld = dr_dt
     
  
  end if


  return
end subroutine Driver_verifyInitDt








