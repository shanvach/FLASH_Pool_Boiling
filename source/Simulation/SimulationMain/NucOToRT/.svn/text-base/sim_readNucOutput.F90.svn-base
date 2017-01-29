!!****if* source/Simulation/SimulationMain/NucOToRT/sim_readNucOutput
!!
!! NAME
!!    sim_readNucOutput
!!
!! SYNOPSIS
!!
!!    sim_readNucOutput(
!!                      integer(IN)    :: part_props,
!!                      integer(IN)    :: maxParticlesPerProc,
!!                      integer(IN)    :: numParticles)
!!
!! DESCRIPTION
!!
!!     Reads in the particles data in a known format from a known file. The number of actual
!!     particles read is returned in the argument numParticles.  This also generates a mapping 
!!     between the FLASH particle property names and the properties from the particle data.
!!
!! ARGUMENTS
!!
!!
!!  part_props : number of properties in the particles datastructure
!!
!!  maxParticlesPerProc : The size of the buffer allocated for particles at compile time. This
!!                        is the largest number of particles that a simulation can have on one
!!                        processor at any time.
!!
!!  numParticles : While coming in it contains the current number of particles mapped to
!!                      this processor. After all the data structure movement, the number
!!                      of local particles might change, and the new value is put back into it
!!
!! NOTES
!!
!!  Particles *will* be moved by the FLASH machinery, and will be alligned to proper blocks by that.
!!  Therefore we can do a simple slicing on the particles and let the movement algorithm take care
!!  of the rest.  This is similar to how the IO unit handles particle reads.
!!
!!  particles : List of particles. It is two dimensional real array, the first dimension
!!              represents each particle's properties, and second dimension is index to
!!              particles.
!!              It is not needed to be an explicit argument in sim_readNucOutput, but will
!!              be allocated in "call io_ptReadParticleData()"
!!***


subroutine sim_readNucOutput(part_props,maxParticlesPerProc,numParticles,iPartFile)


  use IO_data, ONLY : io_chkptFileID

  use Simulation_data, ONLY : sim_nucFileNames,sim_propNames, sim_meshMe,&
       sim_meshNumProcs
  use Logfile_interface, ONLY : Logfile_stamp
  use IO_interface, ONLY : IO_setScalar

  implicit none
#include "constants.h"
#include "Flash_mpi.h"  
  
  integer, intent(IN) :: part_props,maxParticlesPerProc,numParticles,iPartFile

  integer :: fileId, offset, globalNumParticles, localNumParticles, numFilePartProps, i, ierr



  print *,'filename #',iPartFile,': ', sim_nucFileNames(iPartFile)
  

  call HDF5FileOpen(sim_nucFileNames(iPartFile), fileId, 0)

  call getParticleFileInfo(fileId, globalNumParticles, numFilePartProps)
  call IO_setScalar("globalNumParticles", globalNumParticles)

  !DEV: Should probably handle (numFilePartProps > part_props) as an error? - KW

  
  if(numFilePartProps > part_props) then 
     call Logfile_stamp( numFilePartProps, &
          "[sim_readNucOutput] WARNING - Too many particle properties in the file, some will need to be dropped")
  endif
  
  if (.NOT.allocated(sim_propNames))allocate(sim_propNames(numFilePartProps))
  
  call getPropertyNames(fileId, sim_propNames)

  
  !We need a mapping of the read-in particle attributes to what FLASH sees as the particle attributes
  print *, "NumPartProps: ", numFilePartProps, " part_props: ", part_props

  localNumParticles = globalNumParticles/sim_meshNumProcs
  if(mod(globalNumParticles, sim_meshNumProcs) > sim_meshMe) localNumParticles = localNumParticles + 1
  
  offset = sim_meshMe*globalNumParticles/sim_meshNumProcs
  if(mod(globalNumParticles, sim_meshNumProcs) > sim_meshMe) then
     offset = offset + sim_meshMe
  else
     offset = mod(globalNumParticles, sim_meshNumProcs)
  endif


  io_chkptFileID = fileId

  call io_ptReadParticleData()

  call HDF5FileClose(fileId)

  return

end subroutine sim_readNucOutput
