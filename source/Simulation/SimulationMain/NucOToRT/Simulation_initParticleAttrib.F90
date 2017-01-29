!!****if* source/Simulation/SimulationMain/NucOToRT/Simulation_initParticleAttrib
!!
!! NAME
!!
!!  Simulation_initParticleAttrib
!!
!! SYNOPSIS
!!
!!  call Simulation_initParticleAttrib(logical(IN) :: restart)
!!
!! DESCRIPTION
!!
!!  A stub for initializing additional particle attributes if necessary.
!!
!!  This interface is called during initialization after
!!   o  particles positions,
!!   o  particle velocities and other properties needed to advance
!!      particle trajectories,
!!   o  mass properties for active particles,
!!   o  and particle tags
!!  have all been initialized (or possibly read from a checkpoint, if
!!  restart is true).
!!  
!!
!!  By default this does nothing. A typical use would be to initialize
!!  particle properties that were defined in a simulation directory.
!!  It is likely that an implementation only needs to take action
!!  if the restart flag is false.
!!
!!
!! ARGUMENTS
!!
!!  restart - true if restarting from a checkpoint, false otherwise.
!!
!! SEE ALSO
!!
!!  Driver_intFlash
!!  Particles_initPositions
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine Simulation_initParticleAttrib(restart)
  use Grid_interface, ONLY : Grid_mapParticlesToMesh, Grid_getDeltas, Grid_getListOfBlocks,&
       Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr
  use Particles_interface, ONLY : Particles_updateRefinement
  use Particles_data, ONLY:  pt_numLocal, particles, pt_maxPerProc, &
       pt_attributes, &
       pt_meshVar,pt_numAttributes,&
       pt_posInitialized
  use Simulation_data, ONLY:  sim_meshMe, sim_ptMass, theZs, theAs, sim_numAbundanceProps, &
       sim_pProp2A, sim_pProp2Z, sim_iPartFile, sim_ptNumPartFiles, sim_geometry, &
       sim_unkCellWeight, sim_useTrajValues
  use Grid_data, ONLY: gr_geometryOverride, gr_imin

  implicit none
  logical,intent(in) :: restart

  integer       :: i, j, ib, blockID, numBlks
  integer       :: ii, ij, ik
  real,dimension(MDIM) :: deltas
  integer :: blkList(MAXBLOCKS)
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  real :: tmpParticles(pt_numLocal)
  real :: dvol
  integer :: mapMode, pass, npasses

  integer       ::  unk_ind, partProp_ind
  integer       ::  unused_oldLocalNumBlocks

  if (.NOT. restart) then

     call Grid_getListOfBlocks(LEAF, blkList, numBlks)

     if (sim_iPartFile==1) then
        mapMode = 0
     else
        mapMode = 1
     end if

     if (sim_geometry==CYLINDRICAL .AND. gr_geometryOverride .and. &
          gr_imin < 0.0) then
        npasses = 2
     else
        npasses = 1
     end if

     do pass = 1,npasses
#ifdef MASS_PART_PROP
        ! The logic here and elsewhere for counting particles assigned to each cell
        ! assumes that MASS_PART_PROP is 1.0 in all particles!
        !
        ! 
        ! Populating NUMP
        ! 
        partProp_ind = MASS_PART_PROP
        particles(partProp_ind,:) = 1.0
        tmpParticles = particles(partProp_ind,1:pt_numLocal)
        do j = 1,pt_numLocal
           blockID = particles(BLK_PART_PROP,j)
           call Grid_getDeltas(blockID,deltas)
           deltas(NDIM+1:MDIM) = 1.0
           !! Get the volume calculated, to be undone in Particles_mapToMeshOneBlk
           dvol = deltas(IAXIS)*deltas(JAXIS)*deltas(KAXIS)
           particles(partProp_ind,j) = particles(partProp_ind,j) * dvol
        end do
        unk_ind = NUMP_VAR
        if(sim_meshMe==MASTER_PE)print*,'Mapping part prop',partProp_ind,' to UNK var #',unk_ind
        call doMapToMesh(mode=mapMode)
        particles(partProp_ind,1:pt_numLocal) = tmpParticles


        ! 
        ! Populating NUP0, NUMP is exactly the same for now
        ! 
        do j = 1,pt_numLocal
           blockID = particles(BLK_PART_PROP,j)
           call Grid_getDeltas(blockID,deltas)
           deltas(NDIM+1:MDIM) = 1.0
           !! Get the volume calculated, to be undone in Particles_mapToMeshOneBlk
           dvol = deltas(IAXIS)*deltas(JAXIS)*deltas(KAXIS)
           particles(partProp_ind,j) = particles(partProp_ind,j) * dvol
        end do
        unk_ind = NUP0_VAR
        if(sim_meshMe==MASTER_PE)print*,'Mapping part prop',partProp_ind,' to UNK var #',unk_ind
        call doMapToMesh(mode=mapMode)
        particles(partProp_ind,1:pt_numLocal) = tmpParticles

        ! 
        ! Populating PDEN - like NUP0 and NUMP, but scaled differently!
        ! 
        do j = 1,pt_numLocal
           blockID = particles(BLK_PART_PROP,j)
           particles(partProp_ind,j) = particles(partProp_ind,j) * sim_ptMass
        end do
        unk_ind = PDEN_VAR
        if(sim_meshMe==MASTER_PE)print*,'Mapping part prop',partProp_ind,' to UNK var #',unk_ind
        call doMapToMesh(mode=mapMode)
        particles(partProp_ind,1:pt_numLocal) = tmpParticles
#endif

#ifdef NI56_PART_PROP
        if (pass==1) then       !otherwise do it later...
        !
        ! Special Ni56 handling - populate NI56_MSCALAR
        ! 
           partProp_ind = NI56_PART_PROP
           tmpParticles = particles(partProp_ind,1:pt_numLocal)
           particles(partProp_ind,1:pt_numLocal) = particles(partProp_ind,1:pt_numLocal) * sim_pProp2A(partProp_ind) !will undo!
           do j = 1,pt_numLocal
              blockID = particles(BLK_PART_PROP,j)
              call Grid_getDeltas(blockID,deltas)
              deltas(NDIM+1:MDIM) = 1.0
              !! Get the volume calculated, to be undone in Particles_mapToMeshOneBlk
              dvol = deltas(IAXIS)*deltas(JAXIS)*deltas(KAXIS)
              particles(partProp_ind,j) = particles(partProp_ind,j) * dvol
           end do
           unk_ind = NI56_MSCALAR
           if(sim_meshMe==MASTER_PE)print*,'Mapping part prop',partProp_ind,' to UNK var #',unk_ind
           call doMapToMesh(mode=mapMode)
           particles(partProp_ind,1:pt_numLocal) = tmpParticles
        end if
#endif

        if (pass==1) then
           particles(SUMY_PART_PROP,1:pt_numLocal) = 0.0
           particles(SUMZY_PART_PROP,1:pt_numLocal) = 0.0
        end if
        !
        ! Populate all other UNK variables that have one or more particle mappings
        ! in effect (i.e., for which both one or more PARTICLEMAP and particle_attribute_##
        ! statements in Config are in effect).
        ! That is, this loop populates
        !   o  all SPECIES (i.e., variables representing per-element abundances)
        !   o  DENS_VAR, EINT_VAR, TEMP_VAR
        !   o  VEL{X,Y,Z}_VAR
        !   o  RPV{1,2,3}_MSCALAR
        ! 
        do i = 1,pt_numAttributes
           if(pt_meshVar(PT_MAP,i)==PARTICLEMAP_UNK) then
              partProp_ind = pt_attributes(i)
              unk_ind = pt_meshVar(PT_VAR,i)
              if (sim_useTrajValues(unk_ind)) then
                 if (pass==1 .AND. unk_ind == TEMP_VAR) then
                    print*,'Converting TEMP_PART_PROP (for TEMP_VAR) from [10^9 K] to [K]'
                    particles(TEMP_PART_PROP,1:pt_numLocal) = particles(TEMP_PART_PROP,1:pt_numLocal) *1.0e9
                 end if
                 if (pass==1 .AND. &
                      (unk_ind .GE. SPECIES_BEGIN .AND. unk_ind .LE. SPECIES_END)) then
                    particles(SUMY_PART_PROP,1:pt_numLocal) = &
                         particles(SUMY_PART_PROP,1:pt_numLocal) + particles(partProp_ind,1:pt_numLocal)
                    particles(SUMZY_PART_PROP,1:pt_numLocal) = &
                         particles(SUMZY_PART_PROP,1:pt_numLocal) + particles(partProp_ind,1:pt_numLocal) * sim_pProp2Z(partProp_ind)
                    particles(partProp_ind,1:pt_numLocal) = particles(partProp_ind,1:pt_numLocal) * sim_pProp2A(partProp_ind) !will NOT undo!
                 end if
                 tmpParticles = particles(partProp_ind,1:pt_numLocal)
                 do j = 1,pt_numLocal
                    blockID = particles(BLK_PART_PROP,j)
                    call Grid_getDeltas(blockID,deltas)
                    deltas(NDIM+1:MDIM) = 1.0
                    !! Get the volume calculated, to be undone in Particles_mapToMeshOneBlk
                    dvol = deltas(IAXIS)*deltas(JAXIS)*deltas(KAXIS)
                    particles(partProp_ind,j) = particles(partProp_ind,j) * dvol
                 end do
                 if(sim_meshMe==MASTER_PE)print*,'Mapping part prop',partProp_ind,' to UNK var #',unk_ind
                 call doMapToMesh(mode=1)
                 particles(partProp_ind,1:pt_numLocal) = tmpParticles
              else
                 if(sim_meshMe==MASTER_PE)print*,'NOT Mapping part prop',partProp_ind,' to UNK var #',unk_ind,',using chkpt value'
              end if
           end if
        end do


        !
        ! Special SUMY handling - populate SUMY_MSCALAR
        ! 
        partProp_ind = SUMY_PART_PROP
        unk_ind = SUMY_MSCALAR
        if (sim_useTrajValues(unk_ind)) then
           tmpParticles = particles(partProp_ind,1:pt_numLocal)
           do j = 1,pt_numLocal
              blockID = particles(BLK_PART_PROP,j)
              call Grid_getDeltas(blockID,deltas)
              deltas(NDIM+1:MDIM) = 1.0
              !! Get the volume calculated, to be undone in Particles_mapToMeshOneBlk
              dvol = deltas(IAXIS)*deltas(JAXIS)*deltas(KAXIS)
              particles(partProp_ind,j) = particles(partProp_ind,j) * dvol
           end do
           if(sim_meshMe==MASTER_PE)print*,'Mapping part prop',partProp_ind,' to UNK var #',unk_ind
           call doMapToMesh(mode=mapMode)
           particles(partProp_ind,1:pt_numLocal) = tmpParticles
        else
           if(sim_meshMe==MASTER_PE)print*,'NOT Mapping part prop',partProp_ind,' to UNK var #',unk_ind,'(SUMY_MSCALAR)'
        end if

        !
        ! Special SUMZY handling - populate YE_MSCALAR
        ! 
        partProp_ind = SUMZY_PART_PROP
        unk_ind = YE_MSCALAR
        if (sim_useTrajValues(unk_ind)) then
           tmpParticles = particles(partProp_ind,1:pt_numLocal)
           do j = 1,pt_numLocal
              blockID = particles(BLK_PART_PROP,j)
              call Grid_getDeltas(blockID,deltas)
              deltas(NDIM+1:MDIM) = 1.0
              !! Get the volume calculated, to be undone in Particles_mapToMeshOneBlk
              dvol = deltas(IAXIS)*deltas(JAXIS)*deltas(KAXIS)
              particles(partProp_ind,j) = particles(partProp_ind,j) * dvol
           end do
           if(sim_meshMe==MASTER_PE)print*,'Mapping part prop',partProp_ind,' to UNK var #',unk_ind
           call doMapToMesh(mode=mapMode)
           particles(partProp_ind,1:pt_numLocal) = tmpParticles
        else
           if(sim_meshMe==MASTER_PE)print*,'NOT Mapping part prop',partProp_ind,' to UNK var #',unk_ind,'(YE_MSCALAR)'
        end if

#ifdef NI56_PART_PROP
        if (pass==2) then       !deferred to here for pass 2...
        !
        ! Special Ni56 handling - populate NI56_MSCALAR
        ! 
           partProp_ind = NI56_PART_PROP
           tmpParticles = particles(partProp_ind,1:pt_numLocal)
           do j = 1,pt_numLocal
              blockID = particles(BLK_PART_PROP,j)
              call Grid_getDeltas(blockID,deltas)
              deltas(NDIM+1:MDIM) = 1.0
              !! Get the volume calculated, to be undone in Particles_mapToMeshOneBlk
              dvol = deltas(IAXIS)*deltas(JAXIS)*deltas(KAXIS)
              particles(partProp_ind,j) = particles(partProp_ind,j) * dvol
           end do
           unk_ind = NI56_MSCALAR
           if(sim_meshMe==MASTER_PE)print*,'Mapping part prop',partProp_ind,' to UNK var #',unk_ind
           call doMapToMesh(mode=mapMode)
           particles(partProp_ind,1:pt_numLocal) = tmpParticles
        end if
#endif

        if (npasses>1) then
           do j = 1,pt_numLocal
              particles(POSX_PART_PROP,j) = - particles(POSX_PART_PROP,j)
#ifdef POSZ_PART_PROP
              particles(POSZ_PART_PROP,j) = - particles(POSZ_PART_PROP,j)
#endif
              particles(VELX_PART_PROP,j) = - particles(VELX_PART_PROP,j)
              particles(VELZ_PART_PROP,j) = - particles(VELZ_PART_PROP,j)
           end do
           if (pass==1) then
              call Particles_updateRefinement(unused_oldLocalNumBlocks)
              print*,'Simulation_initParticleAttrib: 2nd round follows, local number of particles',pt_numLocal,' on',sim_meshMe
           else if (pass==2 .AND. sim_iPartFile==sim_ptNumPartFiles) then
              call Particles_updateRefinement(unused_oldLocalNumBlocks)
              print*,'Simulation_initParticleAttrib: 2nd round done, particles moved back, local number of particles',&
                   pt_numLocal,' on',sim_meshMe
           end if
        end if
     end do
  end if
contains
  subroutine doMapToMesh(mode)
    use physicaldata,ONLY: unk

    integer,intent(IN) :: mode
    integer :: j, mapMode

    if (pass==1 .AND. sim_unkCellWeight(unk_ind)==0.0) then
       mapMode = mode
    else
       mapMode = 1
    end if
    print*
#ifdef DEBUG_VERBOSE
998 format ('Before: maxval(UNK(',i3,',:,:,1,1)=',1PG23.15)
999 format ('After : maxval(UNK(',i3,',:,:,1,1)=',1PG23.15)
    print 998,unk_ind,maxval(UNK(unk_ind,:,:,1,1))
#endif
    call Grid_mapParticlesToMesh(particles,NPART_PROPS,pt_numLocal,pt_maxPerProc,&
         partProp_ind,unk_ind, mode=mapMode)
#ifdef DEBUG_VERBOSE
    print 999,unk_ind,maxval(UNK(unk_ind,:,:,1,1))
#endif

  end subroutine doMapToMesh
end subroutine Simulation_initParticleAttrib
