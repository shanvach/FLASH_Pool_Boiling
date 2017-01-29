!!****if* source/Particles/ParticlesMain/passive/PredictorCorrectorOld/pt_advancePassive
!!
!! NAME
!!
!!  pt_advancePassive
!!
!! SYNOPSIS
!!
!!  pt_advancePassive(real(in) :: dtOld,
!!                    real(in) :: dtNew,
!!                    real(inout):: particles(:,p_count),
!!                    integer(in):: p_count))
!!
!! DESCRIPTION
!!
!!  Time advancement routine for the passive particle module.
!!  OLD OBSOLETE "PredictorCorrector LOGIC.
!!
!!
!! Problem:
!!  This is only a second-order method when dtNew == dtOld (exactly).
!!
!! ARGUMENTS
!!
!!   dtOld -- previous time interval
!!   dtNew -- current time interval
!!   particles -- particles to advance
!!   p_count  -- the number of particles in the list to advance
!!  
!! PARAMETERS
!!
!!    pt_dtChangeTolerance REAL [0.4] Do Euler step if change in time step is greater than this
!!                                    percentage.  Set to 0 to always do Euler, set to a huge
!!                                    number to always use estimated midpoint velocities.
!!
!! SIDE EFFECTS
!!
!!  Updates the POS{X,Y,Z}, VEL{X,Y,Z}, POSPRED{X,Y,Z}, and VELPRED{X,Y,Z} properties of
!!  particles in the particles structure.
!!
!! NOTES
!!
!!   At the first call, it is assumed that particle initialization
!!   has filled in initial velocity components properly.
!!
!!   This is a slightly corrected version of the scheme that used to be called
!!   "PredictorCorrector" before the FLASH3.0 release.
!!
!!***

!===============================================================================
#ifdef DEBUG_ALL
#define DEBUG_PARTICLES
#endif

subroutine pt_advancePassive (dtOld,dtNew,particles,p_count)
    
  use Particles_data, ONLY: useParticles, pt_typeInfo,&
       pt_posAttrib, pt_velNumAttrib,pt_velAttrib,&
       pt_posPredAttrib, pt_velPredAttrib,&
       pt_dtChangeTolerance, pt_restart

  use Grid_interface, ONLY : Grid_mapMeshToParticles
  
  implicit none

#include "constants.h"  
#include "Flash.h"
#include "Particles.h"
  
  real, INTENT(in)  :: dtOld, dtNew
  integer, INTENT(in) :: p_count
  real,dimension(NPART_PROPS,p_count),intent(INOUT) :: particles

  real          :: timeDiff
  integer       :: i,nstep
  logical, save :: firstCall = .true.  ! initialized on compilation
  logical       :: doCorrector=.false., sameTimeStep
  integer       :: mapType
!! ----------------------------------------------------------------------------

  ! Don't do anything if runtime parameter isn't set
  if (.not.useParticles ) return
  mapType=pt_typeInfo(PART_MAPMETHOD,PASSIVE_PART_TYPE)
  
  ! Testing for predictor step, rather bloated because LBR is bad at logical calls
  sameTimeStep = .false.
  
  timeDiff = abs(dtNew - dtOld)
  if (timeDiff .LT. (pt_dtChangeTolerance * dtOld))  sameTimeStep = .true.
 
  doCorrector = sameTimeStep .AND. ((.NOT.firstCall) .OR. pt_restart)
   
  ! For the first timestep, or if timesteps have changed, we do an Euler step.
  if (.NOT. doCorrector) then
#ifdef DEBUG_PARTICLES
     print*,'Euler only at this time!'
#endif
     ! Map the updated gas velocity field onto the current particle positions to
     ! obtain the updated particle velocities.
     
     call Grid_mapMeshToParticles(particles,&
          NPART_PROPS, BLK_PART_PROP,p_count,&
          pt_posAttrib,pt_velNumAttrib,pt_velAttrib,mapType)

     ! Update the particle positions.
     do i = 1, p_count
        particles(POSX_PART_PROP,i) = particles(POSX_PART_PROP,i) + &
             dtNew * particles(VELX_PART_PROP,i)
        particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + &
             dtNew * particles(VELY_PART_PROP,i)
        particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + &
             dtNew * particles(VELZ_PART_PROP,i)
     enddo
     

         ! end of know-nothing Euler step  
  else   ! time steps are essentially the same, previous values exist
     
     !-------------------------------------------------------------------------------
     ! For subsequent timesteps, we have predicted positions and velocities, so do
     ! an "estimated midpoint" step.
     
     
     ! Map the updated gas velocity field onto the predicted positions to obtain the
     !  predicted particle velocities.
     
     call Grid_mapMeshToParticles(particles,&
          NPART_PROPS, BLK_PART_PROP,p_count,&
          pt_posPredAttrib,pt_velNumAttrib,pt_velAttrib,mapType)

     
     ! Update the particle positions using the predicted particle velocities
     !  (obtained from predicted positions)
     
     do i = 1, p_count
        particles(POSX_PART_PROP,i) =  particles(POSX_PART_PROP,i) + &
             dtNew * 0.5*(particles(VELX_PART_PROP,i) + &
             particles(VELPREDX_PART_PROP,i))
        particles(POSY_PART_PROP,i) = particles(POSY_PART_PROP,i) + &
             dtNew * 0.5*(particles(VELY_PART_PROP,i) + &
             particles(VELPREDY_PART_PROP,i))
        particles(POSZ_PART_PROP,i) = particles(POSZ_PART_PROP,i) + &
             dtNew * 0.5*(particles(VELZ_PART_PROP,i) + &
             particles(VELPREDZ_PART_PROP,i))
     enddo
     ! end of correction section

     ! Map the updated gas velocity field onto the current particle positions to
     ! obtain the updated particle velocities.
     
     call Grid_mapMeshToParticles(particles,&
          NPART_PROPS, BLK_PART_PROP,p_count,&
          pt_posAttrib,pt_velNumAttrib,pt_velAttrib,mapType)

  endif  ! end of split between doing correction and doing Euler timestep

    

  
  ! update initialization
  if (firstCall) then
     firstCall = .false.
  endif
  

  return
  
end subroutine pt_advancePassive


