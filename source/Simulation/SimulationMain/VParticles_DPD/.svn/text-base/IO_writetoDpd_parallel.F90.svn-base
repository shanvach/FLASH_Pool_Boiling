

#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

subroutine IO_writeDpd_parallel(particles,p_count,firstfileflag)
  
  use HDF5 ! This module contains all necessary modules 

  
  use Driver_data, ONLY: dr_globalMe, dr_nbegin, &
       dr_nend, dr_dt, dr_wallClockTimeLimit, &
       dr_tmax, dr_simTime, dr_redshift, &
       dr_nstep, dr_dtOld, dr_dtNew, dr_nbegin, dr_restart
  use Simulation_data, ONLY : domainsize, pt_NumPart
  use Logfile_interface, ONLY : Logfile_stamp
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Particles_interface, ONLY : Particles_getLocalNum, &
       Particles_updateAttributes, Particles_manageLost, Particles_sinkSyncWithParticles
  use Particles_data, ONLY: pt_maxPerProc
  use Grid_data,ONLY : gr_globalDomain, gr_meshNumProcs 
  use Grid_interface, ONLY:Grid_sortParticles
  use IO_data, ONLY : io_meshComm, io_meshMe, io_meshNumProcs,io_globalMe
  use io_ptInterface, ONLY : io_ptWriteParticleData
  use IOParticles_data, ONLY : io_particleFileNumber, &
       io_particleFileIntervalTime, io_nextParticleFileTime, &
       io_particleFileIntervalStep, io_nextParticleFileStep, &
       io_dumpParticleFileExist, useParticles, &
       io_particleFileIntervalZ, io_nextParticleFileZ
#ifdef USE_IO_C_INTERFACE
  use iso_c_binding, ONLY : c_loc
  use io_c_interface, ONLY : io_create_dataset, io_xfer_cont_slab
#else
#define c_loc(x) x
#endif

  implicit none
  

!#include "Flash.h"
!#include "constants.h" 
#include "Flash_mpi.h"  
#include "Particles.h"
  !----------- Arguments list -------------------
  integer,INTENT(INOUT):: p_count
  real,dimension(NPART_PROPS,p_count),INTENT(INOUT)::particles
  logical,optional, INTENT(INOUT) :: firstfileflag
  !-----------------------------------

  ! - Local arguments
  character(len=6) :: index_count,index_proc
  logical :: forceParticleFile 


  integer :: i , ierr

  character (len=MAX_STRING_LENGTH) :: filename
  integer :: pptFileID

 ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(2,2) :: strBuff
  integer           :: particleOffset, localNumParticles
  integer           :: globalNumParticles
  character (len=OUTPUT_PROP_LENGTH) :: partAttributeLabels(NPART_PROPS)

  real          :: lsimTime   !local sim time variable
  real, save    :: oldsimTime = 0.0  !initial value does not matter
  integer       :: lsimGen    !local sim generation variable
  integer, save :: oldsimGen = -1 !different from any valid generation, to trigger update at startup.

  logical, save :: firstCall = .true.
  logical  :: restart, particlesToCheckpoint=.true.
  

  globalNumParticles = pt_NumPart
  localNumParticles = p_count
 
  call int2char(dr_nstep,index_count)
  call int2char(dr_globalMe,index_proc)
  forceParticleFile = .false.
  
  call Timers_start("write_DPD_Particles")
  
  ! open the file
  call io_initFile(io_particleFileNumber, pptFileID, filename, "_part_", forceParticleFile)
  
  
  if (io_globalMe == MASTER_PE) then
     write (strBuff(1,1), "(A)") "type"
     write (strBuff(1,2), "(A)") "particles"
     write (strBuff(2,1), "(A)") "name"
     write (strBuff(2,2), "(A)") trim(filename)
     call Logfile_stamp( strBuff, 2, 2, "[IO_writeParticles] open")
     write(*,*)'This is the particles files name ', trim(filename)
  end if
  
  
  if(oldsimTime /= lsimTime .OR. oldsimGen /= lsimGen) then
     call Particles_updateAttributes()
     if (io_globalMe == MASTER_PE) then
        call Logfile_stamp( 'done called Particles_updateAttributes()', "[IO_writeParticles]")
     end if
     oldsimTime=lsimTime
     oldsimGen =lsimGen
  end if

  ! pull in sink particles
  call Particles_sinkSyncWithParticles(sink_to_part=.true.)


  !map the particle property names from an int to a string
  do i=1, NPART_PROPS
     call Simulation_mapIntToStr(i, partAttributeLabels(i),MAPBLOCK_PART)
  end do

  call Particles_manageLost(PART_EXPAND)

  call io_getParticleOffset( localNumParticles, globalNumParticles, particleOffset)  


  if(globalNumParticles <= 0 .and. io_globalMe == MASTER_PE) then
     print *, "WARNING: globalNumParticles = 0!!!"
  end if
  
  call io_ptWriteParticleData( pptFileID, globalNumParticles, &
       localNumParticles, particleOffset, partAttributeLabels, particlesToCheckpoint)
  
  call io_closeFile( pptFileID)
  
  !increment the particle plotfile number
  !only increment if we are writing a single particle dataset to a file
  io_particleFileNumber = io_particleFileNumber + 1
  
  call Timers_stop("write_DPD_Particles")
  
end subroutine IO_writeDpd_parallel

