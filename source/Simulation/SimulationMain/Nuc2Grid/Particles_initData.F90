!!****if* source/Simulation/SimulationMain/Nuc2Grid/Particles_initData
!!
!! NAME
!!
!!  Particles_initData
!!
!! SYNOPSIS
!!
!!  call Particles_initData(logical(IN) :: restart)
!!
!! DESCRIPTION
!!
!!  Initialization of particle attributes after they have been
!!  dropped into their positions.
!!
!!  This interface is called during Flash initialization after
!!   o  particles positions,
!!   o  mass properties for active particles,  (NOTE: it is true in
!!                                              most cases, but not a 
!!                                              requirement)
!!   o  and particle tags
!!  have been initialized (or possibly read from a checkpoint, if
!!  restart is true).
!!
!!  This is where velocity components are initialized.
!!  For passive particles, this is done by mapping from
!!  the fluid velocity field. The sequence of initialization
!!  is such that fluid fields are guaranteed to have been
!!  initialized to their final starting values, and the
!!  configuration of the grid has been finalized including
!!  creation of blocks up to the maximum refinement level,
!!  when the Particles_initData interface is called.
!!
!!  The default implementation in the Particles main subunit
!!  calls Simulation_initParticleAttrib before returning, so
!!  if custom initialization of additional particle properties
!!  is needed, that's the place to implement it.
!!
!! ARGUMENTS
!!
!!  restart - true if restarting from a checkpoint, false otherwise.
!!
!!
!! SEE ALSO
!!
!!  Driver_intFlash
!!  Simulation_initParticleAttrib
!!
!!***

subroutine Particles_initData(restart, partPosInitialized)
  use Simulation_interface,ONLY : Simulation_initParticleAttrib
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface,ONLY : Grid_mapMeshToParticles, Grid_sortParticles
  use Particles_data,ONLY : particles, pt_numLocal, pt_maxPerProc, &
       pt_posAttrib, pt_velNumAttrib, pt_velAttrib, pt_velInitialized,&
       useParticles, pt_typeInfo,pt_numLocal,&
       pt_posInitialized
  use Simulation_data, ONLY: sim_ptNumPartFiles, sim_iPartFile
  use Particles_interface, ONLY : Particles_initPositions, Particles_finalize
  use pt_interface, ONLY : pt_updateTypeDS
  implicit none
#include "Flash.h"
#include "constants.h"
#include "Particles.h"
  logical, intent(IN) :: restart
  logical, intent(INOUT) :: partPosInitialized
  integer :: part_props=NPART_PROPS
  logical :: updateRefine
  integer :: i,p_begin,p_count,p_end
  integer :: iPartFile
  integer, dimension(MAXBLOCKS,NPART_TYPES) :: particlesPerBlk

  ! Return immediately if useParticles is false.
  if (.NOT. useParticles) then
     return
  end if

  updateRefine = .false.
  !We need a work around for pure active particle simulations.
  !NONEXISTENT is known to be -1, which does not conflict with 
  !real particle types.
  !Note the undefine at the end of this subroutine.


  if (.NOT. restart) then
     do iPartFile = 1,sim_ptNumPartFiles
        if (iPartFile > 1) then
           call Particles_finalize()
           pt_posInitialized = .FALSE.
           pt_velInitialized = .FALSE.
           partPosInitialized = .FALSE.
        end if
        sim_iPartFile = iPartFile
        call Particles_initPositions(partPosInitialized,updateRefine)
        pt_posInitialized = partPosInitialized
#if NPART_TYPES > 1
        call Grid_sortParticles(particles,NPART_PROPS,pt_numLocal,NPART_TYPES, &
             pt_maxPerProc,particlesPerBlk,BLK_PART_PROP,TYPE_PART_PROP)
#else
        call Grid_sortParticles(particles,NPART_PROPS,pt_numLocal,NPART_TYPES, &
             pt_maxPerProc,particlesPerBlk,BLK_PART_PROP)
#endif
        ! Now update the pt_typeInfo data structure
        call pt_updateTypeDS(particlesPerBlk)
        if(.not.partPosInitialized)&
             call Driver_abortFlash("initialization of Particles positions failed")
        if (.NOT. pt_velInitialized) then
           do i = 1,NPART_TYPES
              if((pt_typeInfo(PART_ADVMETHOD,i)==RUNGEKUTTA).or. &
                   (pt_typeInfo(PART_ADVMETHOD,i)==MIDPOINT).or.  &
                   (pt_typeInfo(PART_ADVMETHOD,i)==ESTI).or.&
                   (pt_typeInfo(PART_ADVMETHOD,i)==EULER_MAS))then
                 p_begin=pt_typeInfo(PART_TYPE_BEGIN,i)
                 p_count=pt_typeInfo(PART_LOCAL,i)
                 p_end=p_begin+p_count-1
                 call Grid_mapMeshToParticles(particles(:,p_begin:p_end),&
                      part_props,BLK_PART_PROP, p_count,&
                      pt_posAttrib,pt_velNumAttrib,pt_velAttrib,&
                      pt_typeInfo(PART_MAPMETHOD,i))
              end if
           end do
        end if
        pt_velInitialized = .TRUE.
        
        ! Call Simulation_initParticleAttrib, that's where we do a lot of work.
        call Simulation_initParticleAttrib(restart)
     end do

  else

     pt_velInitialized = .TRUE.
     ! Call Simulation_initParticleAttrib, normally just a stub, to allow setups some customization
     call Simulation_initParticleAttrib(restart)
  end if

end subroutine Particles_initData
