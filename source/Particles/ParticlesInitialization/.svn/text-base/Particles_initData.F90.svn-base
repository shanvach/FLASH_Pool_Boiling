!!****if* source/Particles/ParticlesInitialization/Particles_initData
!!
!! NAME
!!
!!  Particles_initData
!!
!! SYNOPSIS
!!
!!  call Particles_initData(logical(IN)    :: restart,
!!                          logical(INOUT) :: partPosInitialized)
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
!! DEV: The above on order of calls may not be true any more! Check and update text! - KW
!!
!!  This is where velocity components are initialized.
!! DEV: This may also be where Particles posns are initialized (indirectly)!
!!  For passive particles, this is done by mapping from
!!  the fluid velocity field. The sequence of initialization
!!  is such that fluid fields are guaranteed to have been
!!  initialized to their final starting values, and the
!!  configuration of the grid has been finalized including
!!  creation of blocks up to the maximum refinement level,
!!  when the Particles_initData interface is called.
!!
!!  The default implementation in the ParticlesInitialization
!!  subunit calls Simulation_initParticleAttrib before returning,
!!  so if custom initialization of additional particle properties
!!  is needed, that's the place to implement it.
!!
!! ARGUMENTS
!!
!!  restart - true if restarting from a checkpoint, false otherwise.
!!
!! SIDE EFFECTS
!!
!!  Updates particles data that is private to the Particles unit.
!!
!!  Normally sets an internal flag (pt_velInitialized) that keeps track of
!!  whether initialization of particle velocities is complete.
!!
!! DEV: may set internal flag pt_posInitialized (indirectly, by calling Particles_initPositions).
!! DEV: Updates the pt_typeInfo data structure (by calling pt_updateTypeDS).

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
       useParticles, pt_meshMe, pt_meshNumProcs,pt_typeInfo,pt_numLocal,&
       pt_indexCount, pt_indexList
  use Particles_interface, ONLY : Particles_initPositions
  use pt_interface, ONLY : pt_updateTypeDS
  implicit none
#include "Flash.h"
#include "constants.h"
#include "Particles.h"
  logical, intent(IN) :: restart
  logical, intent(INOUT) :: partPosInitialized
  integer :: part_props=NPART_PROPS
  logical :: updateRefine, needVel, coords_in_blk
  integer :: i,p_begin,p_count,p_end
  integer, dimension(MAXBLOCKS,NPART_TYPES) :: particlesPerBlk

  ! Return immediately if useParticles is false.
  if (.NOT. useParticles) then
     return
  end if

  updateRefine = .false.


  if (.NOT. restart) then
     call Particles_initPositions(partPosInitialized,updateRefine)
#ifdef TYPE_PART_PROP
  call Grid_sortParticles(particles,NPART_PROPS,pt_numLocal,NPART_TYPES, &
       pt_maxPerProc,particlesPerBlk,BLK_PART_PROP, TYPE_PART_PROP)
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
           needVel=(pt_typeInfo(PART_ADVMETHOD,i)==RUNGEKUTTA)
           needVel=(pt_typeInfo(PART_ADVMETHOD,i)==ESTI).or.needVel
           needVel=(pt_typeInfo(PART_ADVMETHOD,i)==EULER_TRA).or.needVel
           needVel=(pt_typeInfo(PART_ADVMETHOD,i)==MIDPOINT).or.needVel
           if(needVel) then
              p_begin=pt_typeInfo(PART_TYPE_BEGIN,i)
              p_count=pt_typeInfo(PART_LOCAL,i)
              p_end=p_begin+p_count-1
              call Grid_mapMeshToParticles(particles(:,p_begin:p_end),&
                   part_props,BLK_PART_PROP,p_count,&
                   pt_posAttrib,pt_velNumAttrib,pt_velAttrib,&
                   pt_typeInfo(PART_MAPMETHOD,i))
           end if
        end do
     end if
  end if
  pt_velInitialized = .TRUE.
#ifndef FIXEDBLOCKSIZE
  coords_in_blk=.true.
  if(restart)call Grid_moveParticles(particles,NPART_PROPS,pt_maxPerProc,&
       pt_numLocal,pt_indexList, pt_indexCount, coords_in_blk)
#endif
  ! Call Simulation_initParticleAttrib, normally just a stub, to allow setups some customization
  call Simulation_initParticleAttrib(restart)


end subroutine Particles_initData
