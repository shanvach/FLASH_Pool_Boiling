!!****if* source/Particles/ParticlesMain/passive/PredictorCorrectorOld/pt_preparePassive
!!
!! NAME
!!
!!  pt_preparePassive
!!
!! SYNOPSIS
!!
!!  call pt_preparePassive(real(in)   :: dtOld,
!!                         real(in)   :: dtNew,
!!                         real(inout):: particles(:,p_count),
!!                         integer(in):: p_count))
!!
!! DESCRIPTION
!!
!!  Time advancement routine for the passive particle module.
!!  This portion cleans up and prepares for the next time step.
!!  
!!
!! ARGUMENTS
!!
!!   dtOld -- previous time interval
!!   dtNew -- current time interval
!!   particles -- particles to advance
!!   p_count  -- the number of particles in the list to advance
!!  
!!
!! SIDE EFFECTS
!!
!!  Updates the POSPRED{X,Y,Z} and VELPRED{X,Y,Z} properties of
!!  particles in the particles structure.
!!
!!***

!===============================================================================
#ifdef DEBUG_ALL
#define DEBUG_PARTICLES
#endif

subroutine pt_preparePassive (dtOld,dtNew,particles,p_count)
  
  use Particles_data, ONLY: useParticles, pt_typeInfo,&
  pt_velNumAttrib,pt_posPredAttrib, pt_velPredAttrib
       
       use Grid_interface, ONLY : Grid_mapMeshToParticles
  
  implicit none
  
#include "constants.h"
#include "Flash.h"
#include "Particles.h"
  
  real, INTENT(in)  :: dtOld, dtNew
  
  integer, INTENT(in) :: p_count
  real,dimension(NPART_PROPS,p_count),intent(inout) :: particles
  
  integer :: mapType,i
  
  mapType=pt_typeInfo(PART_MAPMETHOD,PASSIVE_PART_TYPE)
  
  ! Compute the predicted midpoint positions.  These are staggered forward 1/2 step.
  !print *, 'about to compute estimated positions'
  do i = 1, p_count
     particles(POSPREDX_PART_PROP,i) =  particles(POSX_PART_PROP,i) + &
          0.5*dtNew * particles(VELX_PART_PROP,i)
     particles(POSPREDY_PART_PROP,i) =  particles(POSY_PART_PROP,i) + &
          0.5*dtNew * particles(VELY_PART_PROP,i)
     particles(POSPREDZ_PART_PROP,i) =  particles(POSZ_PART_PROP,i) + &
          0.5*dtNew * particles(VELZ_PART_PROP,i)
  enddo
  


  ! Map the current gas velocity field onto the estimated positions to get the
  ! predicted velocities.
  
  call Grid_mapMeshToParticles(particles,NPART_PROPS, BLK_PART_PROP,p_count,&
       pt_posPredAttrib,pt_velNumAttrib,pt_velPredAttrib,mapType)
  
  
  return
  
end subroutine pt_preparePassive

