!!****if* source/Simulation/SimulationMain/Nuc2Grid/io_ptReadParticleData
!!
!! NAME
!!
!! io_ptReadParticleData
!!
!!
!! SYNOPSIS
!!
!! io_ptReadParticleData()
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
!!
!!
!! NOTES
!!
!!
!!***


subroutine io_ptReadParticleData()

  use IO_data, ONLY : io_outputSplitNum, &
       io_chkptFileID, io_comm, io_splitNumBlks, io_splitParts, &
       io_meshMe, io_meshNumProcs
  use IOParticles_data, ONLY: io_ptMaxReadPerProc
  use Simulation_interface, ONLY : Simulation_mapIntToStr, Simulation_mapStrToInt
  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stamp
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Particles_interface, ONLY : Particles_putLocalNum
  use Grid_interface, ONLY : Grid_getLocalNumBlks

  use Grid_data, ONLY : gr_globalNumBlocks, gr_globalOffset
  use tree, ONLY : nodetype
  use Simulation_data, ONLY: theZs, theAs, sim_numAbundanceProps, theSpeciesNames, &
       sim_pProp2A, sim_pProp2Z

  use Particles_data, ONLY : particles, pt_numLocal, pt_maxPerProc, pt_logLevel
  use IO_interface, ONLY : IO_getScalar, IO_setScalar


  implicit none

#include "constants.h"
#include "Flash_mpi.h"
#include "Flash.h"
#include "Particles.h"



  integer           :: partsToLeft(0:MAX_PROCS-1)
  integer :: ierr, i, particleOffset
  integer, save :: globalNumParticles !done for IBM compilers
  integer :: localNumBlocks, blkid, lb, startIndex, endIndex

  integer :: reLocalNumParticles  !restart local number of particles, could be different from 
  !what was written out
  integer,allocatable :: particlesPerBlk(:)

  logical :: useParticles, outside

  real,dimension(2,MDIM) :: bndBox !for debugging only ...

  !for file splitting
  integer :: localOffset, splitOffset, partSplitOffset, localPartOffset 

  !for property-by-property read-in.
  integer :: fileNumPartProps, j
  character(len=24), allocatable :: filePropNames(:) 
!!$  integer :: fileToCurrentMap(NPART_PROPS)
  integer,allocatable :: fileToCurrentMap(:)
  integer :: propIndex, pLeft
  character(len=24) :: propString
  logical :: allParticlePropsSame

  ! for selecting particles -Min
  logical :: doRandomSelect
  real :: particlesFraction
  real, allocatable, dimension(:,:) :: particles2
  integer :: jpt, particlesMod, numSelect, last, ith
  real, dimension(10) :: randNum


  !need to get previous runtime parameter storing max particles per proc
  call RuntimeParameters_get("pt_maxPerProc", pt_maxPerProc)
  print*,'io_ptReadParticleData: pt_maxPerProc is', pt_maxPerProc

  call RuntimeParameters_get("useParticles", useParticles)

  if(.not. useParticles) then
     return
  end if


  !allocate particles data structure
  if (.NOT.allocated(particles)) then
     allocate (particles(NPART_PROPS,pt_maxPerProc), stat=ierr)
     if (ierr /= 0) &
          call Driver_abortFlash("io_ptReadParticleData:  could not allocate particle array")
  end if

  !particles must be initialized or the entire particles algorithm will fail
  particles = NONEXISTENT

  call IO_getScalar("globalNumParticles", globalNumParticles)
  if (io_outputSplitNum /= 1) then
     call IO_getScalar("splitNumParticles", io_splitParts)
     call IO_getScalar("splitNumBlocks", io_splitNumBlks)
  else
     io_splitParts = globalNumParticles
     io_splitNumBlks = gr_globalNumBlocks
  endif

  call Grid_getLocalNumBlks(localNumBlocks)

  call MPI_ALLREDUCE(gr_globalOffset, splitOffset, 1, FLASH_INTEGER, &
       MPI_MIN, io_comm, ierr);
  localOffset = gr_globalOffset - splitOffset

  reLocalNumParticles = 0


  if(globalNumParticles > 0 ) then

     allocate(particlesPerBlk(localNumBlocks))

     particlesPerBlk(:) = io_splitParts !complete rubbish...
     particlesPerBlk(:) = 0

     reLocalNumParticles = globalNumParticles/io_meshNumProcs
     if(mod(globalNumParticles, io_meshNumProcs) > io_meshMe) reLocalNumParticles = reLocalNumParticles + 1

     pLeft = reLocalNumParticles
     j = 0
     do while (pLeft > 0 .AND. j .LE. localNumBlocks)
        j = j+1
        i = 1
        do while (pLeft > 0 .AND. i .LE. localNumBlocks)
           if (nodetype(i) == LEAF) then
              particlesPerBlk(i) = particlesPerBlk(i) + 1
              pLeft = pLeft - 1
           end if
           i = i + 1
        end do
     end do


     
     if (io_ptMaxReadPerProc < 0) then
        if (reLocalNumParticles > pt_maxPerProc) then
           call Driver_abortFlash &
                ('[io_ptReadParticleData] ERROR: too many particles on this proc; increase pt_maxPerProc')
        end if
     else
        if (io_ptMaxReadPerProc > pt_maxPerProc) then
           call Driver_abortFlash &
                ('[io_ptReadParticleData]' // &
                'ERROR: too many particles on this proc; increase pt_maxPerProc, or lower io_ptMaxReadPerProc')
        end if
     end if

     !now get the particle offset
     localPartOffset = io_meshMe*(globalNumParticles/io_meshNumProcs)
     if(mod(globalNumParticles, io_meshNumProcs) > io_meshMe) then
        localPartOffset = localPartOffset + io_meshMe
     else
        localPartOffset = localPartOffset + mod(globalNumParticles, io_meshNumProcs)
     endif

     if (io_ptMaxReadPerProc .GE. 0) then
        if (io_ptMaxReadPerProc < reLocalNumParticles) then
           if (pt_logLevel > PT_LOGLEVEL_WARN_DATA) then
              print*,'[io_ptReadParticleData] Number of particles kept',' on PE',io_meshMe,&
                   ' is', io_ptMaxReadPerProc,' out of',reLocalNumParticles
           end if
           reLocalNumParticles = min(io_ptMaxReadPerProc, reLocalNumParticles)
           if (pt_logLevel > PT_LOGLEVEL_WARN_DATA) then
              call Logfile_stamp(reLocalNumParticles,'[io_ptReadParticleData] Number of particles kept')
           end if
        end if
     end if

     !DEV: Changes begin here.  --PR
     !grab the number of particles in the file.

     call io_h5read_num_props(io_chkptFileID, fileNumPartProps)

     allocate(filePropNames(fileNumPartProps))
     allocate(fileToCurrentMap(fileNumPartProps))
     fileToCurrentMap = NONEXISTENT


     sim_numAbundanceProps = fileNumPartProps;
     allocate(theSpeciesNames(fileNumPartProps))

     call io_h5read_globalstrvect(io_chkptFileID, theSpeciesNames, "species names"//char(0), sim_numAbundanceProps)
     do j = 1, sim_numAbundanceProps
        do i=2,len(theSpeciesNames(j))
           if (ichar(theSpeciesNames(j)(i:i)) == 0) theSpeciesNames(j)(i:i) = ' '
        end do
     end do
!!$     print*,'There were',sim_numAbundanceProps,' species names in the file:',theSpeciesNames
     call io_h5read_globalintvect(io_chkptFileID, theAs, "species A"//char(0), sim_numAbundanceProps)
!!$     print*,'There were',sim_numAbundanceProps,' A values for species in the file:',theAs
     call io_h5read_globalintvect(io_chkptFileID, theZs, "species Z"//char(0), sim_numAbundanceProps)
!!$     print*,'There were',sim_numAbundanceProps,' Z values for species in the file:',theZs


     
     do j = 1, sim_numAbundanceProps
        call Simulation_mapStrToInt(theSpeciesNames(j), i, MAPBLOCK_PART)
993     format(' The species (from file) "',A,'", A=',I4,',Z=',I4,' maps to particle property',I4)
        print 993,trim(theSpeciesNames(j)),theAs(j),theZs(j),i

        if (i < 1) then
           call Driver_abortFlash('Invalid particle property for a species!')
        end if
        sim_pProp2A(i) = theAs(j)
        sim_pProp2Z(i) = theZs(j)
     end do

     deallocate(theSpeciesNames)


     call io_h5read_particle_names(io_chkptFileID, filePropNames, fileNumPartProps);

991  format (I3,':"',A,'"')
     do j = 1, fileNumPartProps
        do i=2,len(filePropNames(j))
           if (ichar(filePropNames(j)(i:i)) == 0) filePropNames(j)(i:i) = ' '
        end do
!!        print 991,j,filePropNames(j)
     end do

     !generate mapping
     do i = 1,NPART_PROPS

        !iterate over part
        call Simulation_mapIntToStr(i, propString, MAPBLOCK_PART)
        do j = 1, fileNumPartProps

992        format (I3,':"',A,'"',A2,'"',A,'"')
           if(trim(propString) .eq. trim(filePropNames(j))) then
!!              print 992,j,propstring,'==',filePropNames(j)
              fileToCurrentMap(j) = i
              exit
!!$           else
!!$              print 992,j,propstring,'NE',filePropNames(j)
           end if

        end do


     end do


     !Test whether the particle attributes are identical in the 
     !current simulation & in the checkpoint file.
     allParticlePropsSame = .false.
     if (NPART_PROPS == fileNumPartProps) then
        do i = 1, NPART_PROPS
           if (fileToCurrentMap(i) == i) then
              allParticlePropsSame = .true.
           else
              allParticlePropsSame = .false.
              exit
           end if
        end do
     end if


     !We discovered a strange deadlock on BG/P during the attribute 
     !by attribute particle read.  As such, we now perform the simpler 
     !contiguous data read when possible.  Otherwise we give a very 
     !verbose output during the attribute by attribute read... until we 
     !understand what is happening.
     if (allParticlePropsSame .eqv. .true.) then
        !DEV: May want to keep this for performance reasons...
        if (io_meshMe == MASTER_PE) then
           print *, "[io_ptReadParticleData]: Starting contiguous particle read:"
        end if
        call io_h5read_particles(io_chkptFileID, &
             particles, &
             reLocalNumParticles, &
             NPART_PROPS, &
             localPartOffset)
     else
        if (io_meshMe == MASTER_PE) then
           print *, "[io_ptReadParticleData]: Starting attribute by attribute particle read:"
           print *, "[io_ptReadParticleData]: (io_comm == MPI_COMM_WORLD):", & 
                (io_comm == MPI_COMM_WORLD)
        end if

        do i = 1, fileNumPartProps

           if(fileToCurrentMap(i) .eq. NONEXISTENT ) then
             if (io_meshMe == MASTER_PE) then 
                  print *, "[io_ptReadParticleData]: Skip chkpoint attribute:", &
                       i, "...this is not in the simulation"
              end if
              !force iteration
              cycle
           end if

           if (io_meshMe == MASTER_PE) then
              print *, "[io_ptReadParticleData]: Reading data from chkpoint attribute:", &
                   i, "and storing in simulation attribute:", fileToCurrentMap(i)
           end if
           !print *, "prop map:", i, fileToCurrentMap(i)
           ! print *, fileNumPartProps, NPART_PROPS

#ifdef DEBUG_GRIDPARTICLES
           print*,'PE',io_meshMe,' about to read',reLocalNumParticles,'@',localPartOffset
#endif
           call io_h5read_single_part_prop(io_chkptFileID, &
                particles, &
                reLocalNumParticles, &
                fileToCurrentMap(i), &
                i, &
                fileNumPartProps, &
                NPART_PROPS, &
                localPartOffset)
           call MPI_BARRIER(io_comm, ierr)
           !print *, particles(i,:)
        end do

     end if


     deallocate(fileToCurrentMap)


     ! select a subset of paticles 

     call RuntimeParameters_get("doRandomSelect", doRandomSelect)
     call RuntimeParameters_get("particlesFraction", particlesFraction)
     call RuntimeParameters_get("particlesMod", particlesMod)

     allocate(particles2(NPART_PROPS, reLocalNumParticles))
     
     particles2(:,1:reLocalNumParticles) = particles(:,1:reLocalNumParticles)
     particles = NONEXISTENT

     if(doRandomSelect) then
       numSelect=reLocalNumParticles*particlesFraction
       last=reLocalNumParticles
       call init_random_seed()
       do i=1, numSelect
          call random_number(randNum)
          ith=(last-1)*randNum(10)+1
          particles(:,i)=particles2(:,ith)
          particles2(:,ith)=particles2(:,last)
          last=last-1
       end do
       reLocalNumParticles=numSelect
     else
       jpt  = 0
       do i = 1, reLocalNumParticles
          if (mod(i,particlesMod) == 0) then
            jpt = jpt + 1
            particles(:,jpt) = particles2(:,i)
          end if
       end do

       ! notice new number of selected particles
       reLocalNumParticles = jpt
     end if

     deallocate(particles2)
     ! end of selecting

     if (io_meshMe == MASTER_PE) then
        print *, "[io_ptReadParticleData]: Finished reading particles from file."
     end if


#if (0)
     ! We leave it to later code to
     ! 1) deal with invalid BLK attributes (which may signal a particle that has
     !    left the original domain, and
     ! 2) set BLK attributes to correct values for all particles that we retain
     !    in the process of moving them to the right blocks (and perhaps proc),
     !    see Grid_moveParticlesGlobal call. - KW
     startIndex = 1
     do lb=1, localNumBlocks

        if(particlesPerBlk(lb) > 0) then
           endIndex = startIndex + particlesPerBlk(lb)

           particles(BLK_PART_PROP,startIndex:endIndex-1) = lb
           startIndex = endIndex

        end if
     end do
#endif

     deallocate(filePropNames)

     deallocate(particlesPerBlk)


  end if !if globalNumParticles > 0

  call Particles_putLocalNum(reLocalNumParticles)

#ifdef DEBUG_GRIDPARTICLES
  !check to see if all the particles have a valid BLK_PART_PROP
  do i=1, reLocalNumParticles

     blkid = int(particles(BLK_PART_PROP,i))

     if((blkid < 0) .or. (blkid > localNumBlocks)) then
        call Driver_abortFlash("io_pteadParticleData, blkid out of bounds")
     end if

  end do
#endif

#if(0)
  do i=1, reLocalNumParticles

     blkid = int(particles(BLK_PART_PROP,i))
!!$     if (particles(POSX_PART_PROP,i) .GE. 6.18840E+09  .AND. particles(POSX_PART_PROP,i) .LE. 6.53220E+09 .AND. &
!!$         particles(POSY_PART_PROP,i) .GE. -1.10016E+10 .AND. particles(POSY_PART_PROP,i) .LE. -1.06578E+10) then
     if (particles(POSX_PART_PROP,i) .GE. 5.84460E+09 .AND. particles(POSX_PART_PROP,i) .LE. 6.87600E+09 .AND. &
         particles(POSY_PART_PROP,i) .GE. -1.13454E+10 .AND. particles(POSY_PART_PROP,i) .LE. -1.03140E+10) then
        print*,'********************* This particles is IT! **********************',i,blkid
995     format((1x,'Prop',I3,':',1PG23.16))
        print 995,(j,particles(j,i),j=1,NPART_PROPS)
        print*,'******************************************************************'
     end if


  end do
#endif

  return

end subroutine io_ptReadParticleData


SUBROUTINE init_random_seed()
   INTEGER :: i, n, clock
   INTEGER, DIMENSION(:), ALLOCATABLE :: seed

   CALL RANDOM_SEED(size = n)
   ALLOCATE(seed(n))

   CALL SYSTEM_CLOCK(COUNT=clock)

   seed = clock + 37 * (/ (i - 1, i = 1, n) /)
   CALL RANDOM_SEED(PUT = seed)

   DEALLOCATE(seed)
 END SUBROUTINE

