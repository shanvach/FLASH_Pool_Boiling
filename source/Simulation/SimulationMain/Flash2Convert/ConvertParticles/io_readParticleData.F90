!!****if* source/Simulation/SimulationMain/Flash2Convert/ConvertParticles/io_readParticleData
!!
!! NAME
!!
!! io_readParticleData
!!
!!
!! SYNOPSIS
!!
!! io_readParticleData()
!!                      integer(in) :: numProcs)
!!
!!
!!
!! DESCRIPTION
!!
!!    This routine reads out the particle data in a separate hdf5 file
!!    It calls  ioh5_read_particles
!!    This is being done to make particles easy to debug and since the particles
!!    are not associated with the mesh data it makes since to separate it from
!!    the checkpoint files and plotfiles
!!
!! ARGUMENTS
!!
!!  myPE : current processor number
!!
!!  numProcs : number of processors in the simulation
!!
!!
!! NOTES
!!
!!
!!***


subroutine io_readParticleData(myPE, numProcs)

  use IO_data, ONLY : io_outputSplitNum, &
       io_chkptFileID
  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, RuntimeParameters_getPrev
  use Particles_interface, ONLY : Particles_putLocalNum
  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
    Grid_getListOfBlocks, Grid_sortParticles

  use Grid_data, ONLY : gr_globalNumBlocks, gr_globalOffset

  use Particles_data, ONLY : particles, pt_maxPerProc
  use IO_interface, ONLY : IO_setScalar, IO_getScalar
  use Simulation_interface, ONLY : Simulation_mapStrToInt
  use Grid_interface, ONLY : Grid_sortParticles

  implicit none

#include "constants.h"
#include "mpif.h"
#include "Flash.h"


  integer, intent(in) :: myPE, numProcs

  integer           :: partsToLeft(0:MAX_PROCS-1)
  integer :: ierr, globalNumParticles, i, particleOffset
  integer :: localNumBlocks, blkid, lb, startIndex, endIndex

  integer :: reLocalNumParticles  !restart local number of particles, could be different from 
                                  !what was written out
  integer,allocatable :: particlesPerBlk(:)

  logical :: useParticles, outside

  real,dimension(2,MDIM) :: bndBox !for debugging only ...
  
  !real,allocatable :: tmpParts(:,:)
  integer, allocatable :: localnp(:)
  integer :: key

  character(len=24),allocatable :: intPropNames(:), realPropNames(:)
  integer, allocatable :: intPropIndexMap(:), realPropIndexMap(:)
  character(len=24),allocatable :: propNames(:,:)
  character :: intOrReal
  
  integer :: j,k,numFlash2vars, numIntParts, numRealParts, partx, party, partz
  integer, allocatable :: flash2IndexMap(:), intPropMap(:), realPropMap(:)

  type particle_type
    integer, dimension(2):: intProperty
    real, dimension(NPART_PROPS-2) :: realProperty
  end type particle_type
  type(particle_type), allocatable :: f2_particles(:)

  integer :: blkIdx, blkCount
  integer :: blkList(MAXBLOCKS)

  !information for the file itself
  integer :: fileNumIntProps, fileNumRealProps, localNumParticles
  character (len=24), allocatable :: fileIntPropNames(:), fileRealPropNames(:)
  integer :: partsPerPE(numProcs), localNPBuf(1)
  


  !need to get previous runtime parameter storing max particles per proc
  call RuntimeParameters_get("pt_maxPerProc", pt_maxPerProc)


  call RuntimeParameters_get("useParticles", useParticles)

  if(.not. useParticles) then
     return
  end if

  allocate(particlesPerBlk(MAXBLOCKS))
  !allocate particles data structure
  allocate (particles(NPART_PROPS,pt_maxPerProc), stat=ierr)
  if (ierr /= 0) then
     call Driver_abortFlash("io_readParticleData:  could not allocate particle array")
  endif

  call Grid_getLocalNumBlks(localNumBlocks)

  !PARTICLE MAPPING FROM SEPARATE FILE:

  OPEN(UNIT=10, FILE='flash2parts', STATUS='OLD')
  
  READ(10, *)numFlash2vars, numIntParts, numRealParts
  print('(I4, I4, I4)'), numFlash2vars, numIntParts, numRealParts

  !print *,'checking header info', numFlash2vars, numIntParts, numRealParts 
 
  if(numFlash2vars /= NPART_PROPS) then
    call Driver_abortFlash("[io_readParticleData]: flash2parts total properties do not agree.")
  endif

  allocate(propNames(numFlash2vars, 2))
  allocate(intPropNames(numIntParts))
  allocate(realPropNames(numRealParts))
  !allocate(intPropMap(numIntParts))
  !allocate(realPropMap(numRealParts))

  
  
  j=1
  k=1
  do i = 1, numFlash2vars
    READ(10, *) propNames(i,1), propNames(i,2), intOrReal
    if(intOrReal == 'r') then
       realPropNames(j)=propNames(i,1)
       j=j+1
    else if (intOrReal == 'i') then
      intPropNames(k)=propNames(i,1)
      k=k+1
    else
      call Driver_abortFlash('[io_readParticleData]: Particle property not&
           labeled an int or a real!')
    end if
    print *, propNames(i,1), propNames(i,2), intOrReal
  end do
  CLOSE(10)
  !DEV: should read in on one proc and distribute.

!print *, intPropNames(:), realPropNames(:)

  
  !END FLASH2 TO FLASH3 PARTICLE MAPPING.

  !particles must be initialized or the entire particles algorithm will fail
  particles(:,:) = NONEXISTENT


  call io_get_numparticles(io_chkptFileID, globalNumParticles)

  !call IO_setScalar("globalNumParticles", globalNumParticles)
  call io_setPrevScalarInt("globalNumParticles", globalNumParticles)
  !call IO_getScalar("globalNumParticles", globalNumParticles)
  print *, globalNumParticles

  call Grid_getLocalNumBlks(localNumBlocks)
  
  reLocalNumParticles = 0


  if(globalNumParticles > 0 ) then
     

     !allocate(particlesPerBlk(localNumBlocks))

     

     !TODO: FAKE if needed.
     !DEV: It's needed.
     !**FLASH 2 has no localnp, so have to fake this.**!
     !return an array particlePerBlk holding the number 
     !of particles on each blk on the local proc
     !call io_h5read_localnp(io_chkptFileID, &
     !     particlesPerBlk, &
     !     localNumBlocks, &
     !     gr_globalNumBlocks, &
     !     gr_globalOffset)

     !Faking up local np, a particles move should then get things onto the
     !right processors.  Going with an approximately even distribution.
     localNumParticles = INT(globalNumParticles / numProcs)
     i = mod(globalNumParticles, REAL(numProcs))
     if(myPE < i) localNumParticles = localNumParticles + 1 
     localNPBuf(1) = localNumParticles

     !we will need to figure out what the offset is.
     !DEV: should this be modified for communicator-based IO?
     call MPI_ALLGATHER(localNPBuf, 1, MPI_INTEGER, &
                        partsPerPE, 1, MPI_INTEGER, MPI_COMM_WORLD,ierr)

     print *,partsPerPE
  
     particleOffset = 0
     if(myPE .NE. MASTER_PE) particleOffset = SUM(partsPerPE(1:myPE))
     print *, particleoffset

     
     !DEV:if this is going to scale to large numbers of particles need to just
     !read in localnp of particles (Talk to ZuHone for a test case)
     !split particles somewhat evenly. Must go through a redistribute anyway
     
     !we need the names in the file to intelligently read things into the 
     !particles array, will slicing work here?
     call h5_read_f2_numprops(io_chkptFileID, fileNumIntProps, fileNumRealProps)
     print *, fileNumIntProps, fileNumRealProps
     allocate(fileIntPropNames(fileNumIntProps))
     allocate(fileRealPropNames(fileNumRealProps))
     call sim_h5read_f2PartPropNames(io_chkptFileID, fileIntPropNames, fileRealPropNames)

     allocate(intPropMap(fileNumIntProps))
     allocate(realPropMap(fileNumRealProps))


     !Construct property mapping
     intPropMap(:) = NONEXISTENT
     do i = 1, fileNumIntProps
        do j = 1, numFlash2Vars
           if(fileIntPropNames(i) == propNames(j,1)) then
           print *, fileIntPropNames(i), propNames(j,1)
           call Simulation_mapStrToInt(propNames(j,2),intPropMap(i),&
                MAPBLOCK_PART)
           EXIT
           end if
        end do
     end do

     realPropMap(:) = NONEXISTENT
     do i = 1, fileNumRealProps
        do j = 1, numFlash2Vars
           if(fileRealPropNames(i) == propNames(j,1))then
              print *, fileRealPropNames(i), propNames(j,2)
              call Simulation_mapStrToInt(propNames(j,2),realPropMap(i), MAPBLOCK_PART)
              EXIT
           end if
        end do
     end do

     do i = 1, fileNumIntProps
        print *, fileIntPropNames(i), intPropMap(i)
     end do
 
     do i = 1, fileNumRealProps
        print *,fileRealPropNames(i), realPropMap(i)
     end do
     

     !have to go through the particle tags and find how many parts are on 
     !each block, so read in the whole array, and figure what our localnps are
     
     allocate(localnp(gr_globalNumBlocks))

!    I have to break this up into independent reads. 

     !read in int quantities
     print *, localNumParticles, particleOffset,fileNumIntProps, "****!!!****"
     do i = 1, fileNumIntProps
        if(intPropMap(i) .NE.  NONEXISTENT) then
           print*, fileIntPropNames(i)
           call sim_h5read_f2intProp(io_chkptFileID,&
                                     fileIntPropNames(i), &
                                     localNumParticles, &
                                     particleOffset, &
                                     particles(intPropMap(i),:))
        print *, particles(intPropMap(i),:)
        end if
     end do
     
     do i = 1, fileNumRealProps
        if(realPropMap(i) .NE.  NONEXISTENT) then
           print*, filerealPropNames(i)
           call sim_h5read_f2realProp(io_chkptFileID,&
                                     fileRealPropNames(i), &
                                     localNumParticles, &
                                     particleOffset, &
                                     particles(realPropMap(i),:))
        print *, particles(realPropMap(i),:)
        end if
     end do
     


!!$     allocate(f2_particles(globalNumParticles))
!!$     
!!$     !must read in global particle data for reconstruction.
!!$     call io_h5read_f2_particles(io_chkptFileID,&
!!$       globalNumParticles, &
!!$       globalNumParticles,&
!!$       0, &
!!$       intPropNames,&
!!$       realPropNames,&
!!$       f2_particles)
!print *, 'Particle read OK'
!print *, intPropNames(:), realPropNames(:)
     !convert all particles into a plain double array:
     print *, 'read OK'
     print *, particles
     !print *, f2_particles(1)%intProperty(:)
     
     !Get particle names alone.
     


!!$     do i = 1, globalNumParticles
!!$       do j = 1, numIntParts
!!$     particles(flash2IndexMap(j),i) = REAL(f2_particles(i)%intProperty(j))
!!$       end do
!!$       do j = numIntParts+1, numRealParts+numIntParts
!!$         particles(flash2IndexMap(j),i) = f2_particles(i)%realProperty(flash2IndexMap(j))
!!$       end do
!!$     end do
    
      !print *, particles(1, :)
    
  !   particlesPerBlk(1:localNumBlocks) = localnp(gr_globalOffset:&
!       &gr_globalOffset + localNumBlocks)
  
     deallocate(localnp)
     !deallocate(f2_particles) !Should this be here?
     
     !find the index of a particle's block
     call Simulation_mapStrToInt('blk', blkIdx, MAPBLOCK_PART)
     call Grid_getListOfBlocks(ALL_BLKS, blkList, blkCount)
     print *, blkList
     !now find the newLocalNumParticles
     
     !call gr_ptFindBlock(blkList, blkCount, particles, blkID)
    

     particlesPerBlk(:) = 0
     do lb = 1,blkCount
       do i = 1, globalNumParticles
        if(blkList(lb) == particles(blkIdx, i))then
      particlesPerBlk(lb) =  particlesPerBlk(lb)+1
    end if
       end do
     end do
    
     !print *,particlesPerBlk(:) 
     
     
     do lb=1, localNumBlocks
        reLocalNumParticles = particlesPerBlk(lb) + reLocalNumParticles
     end do

     !print *, reLocalNumParticles 
     if (reLocalNumParticles > pt_maxPerProc) then
        call Driver_abortFlash &
             ('[io_readParticleData] ERROR: too many particles on this proc; increase pt_maxPerProc')
     end if


     !now get the particle offset
     call io_getParticleOffset(myPE, numProcs, reLocalNumParticles, globalNumParticles, particleOffset)


     
     !call io_h5read_particles(io_chkptFileID, &
      !    particles, &
       !   reLocalNumParticles, &
        !  NPART_PROPS, &
         ! particleOffset)

     
     !reset particles BLK_PART_PROP because it could have changed on restart
     startIndex = 1
     do lb=1, localNumBlocks

        if(particlesPerBlk(lb) > 0) then
           endIndex = startIndex + particlesPerBlk(lb)
           
           particles(BLK_PART_PROP,startIndex:endIndex-1) = lb
           startIndex = endIndex
           
        end if
     end do

  end if !if globalNumParticles > 0

  call Particles_putLocalNum(reLocalNumParticles)

  !call gr_packParticles(particles, relocalNumParticles, pt_maxPerProc)
  call Grid_sortParticles(particles, NPART_PROPS, relocalNumParticles, & 
       pt_maxPerProc, particlesPerBlk, .true.)


#ifdef DEBUG_GRIDPARTICLES
  !check to see if all the particles have a valid BLK_PART_PROP
  do i=1, reLocalNumParticles

     blkid = int(particles(BLK_PART_PROP,i))

     if((blkid < 0) .or. (blkid > localNumBlocks)) then
        call Driver_abortFlash("io_readParticleData, blkid out of bounds")
     end if

!!$     !do an expensive check here to see if all particles are within
!!$     call Grid_getBlkBoundbox(int(particles(BLK_PART_PROP,i)), bndBox)
!!$     call gr_particleOutsideBndBox(bndBox, particles(POSX_PART_PROP:POSZ_PART_PROP, i), outside)
!!$     if(outside) then
!!$        print *, "particle outside bndbox ", i
!!$        call Driver_abortFlash("Error: io_readParticleData, particle outside bndBox")
!!$     end if

   end do
#endif

deallocate(propNames)
deallocate(intPropNames)
deallocate(realPropNames)
deallocate(particlesPerBlk)


!call Driver_abortFlash("[io_readParticleData]: diagnostic terminated OK.")
  print *, 'Read Particle Data exiting'
  return

end subroutine io_readParticleData
