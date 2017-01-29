!!****if* source/Simulation/SimulationMain/ShockCyl/Particles_updateAttributes
!!
!! NAME
!!
!!  Particles_updateAttributes
!!
!! SYNOPSIS
!!
!!  Particles_updateAttributes()
!!
!! DESCRIPTION
!!
!!  Particle attribute advancement routines.  Use this routine to 
!!    update user-define particle attributes beyond the usual
!!    particle attributes of position, velocity, block, tag, and mass.
!!
!!  Called indirectly from IO_output.
!!
!!  This simulation-specific implementation has hardwired mappings
!!  and does not use runtime parameters particle_attribute_1, particle_attribute_2,
!!  etc.
!!
!! ARGUMENTS
!!  
!! PARAMETERS
!!
!! NOTES
!!
!!  Added attributes: 
!!    Field variables --  pressure, density, and sf6 mass fraction.
!!    Sim parameters -- nstep, simtime
!!    many more      -- see derivedVariables subroutine called by this implementation.
!!  No filtering.
!!
!! SEE ALSO
!!  sim_derivedVariables
!!***

! DEV this routine needs fixing for when Particles are/are not included in the setup

subroutine Particles_updateAttributes()

  use Particles_data, ONLY:  pt_numLocal, particles, useParticles, pt_maxPerProc,pt_myPe
  use Driver_interface, ONLY : Driver_getNStep, Driver_getSimTime
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_interface, ONLY : Grid_fillGuardCells, Grid_mapMeshToParticles

  implicit none

#include "constants.h"
#include "Flash.h"

  integer       :: i, nstep     
  real          :: simtime              
  integer,parameter :: particleTypes=1

  !------------------------------------------------------------------------------

  ! boot if not using particles
  if (.not. useParticles) return

  !------------------------------------------------------------------------------

  call Driver_getNStep(nstep)
  call Driver_getSimTime(simtime)

  do i = 1, pt_numLocal
     particles(NSTEP_PART_PROP,i) = real(nstep)
     particles(SIMTIME_PART_PROP,i) = simtime
  enddo

  call Timers_start ("Particles_updateAttributes: map from mesh")

  ! Prepare guard cells for mapping mesh variables to particle properties
  call Grid_fillGuardCells(CENTER,ALLDIR,doEos=.true.)
  call Grid_mapMeshToParticles(particles, NPART_PROPS, BLK_PART_PROP,pt_numLocal,particleTypes,pt_maxPerProc,&
       DENS_VAR, DENS_PART_PROP,   &
       POSX_PART_PROP,POSY_PART_PROP,POSZ_PART_PROP)
  call Grid_mapMeshToParticles(particles, NPART_PROPS, BLK_PART_PROP,pt_numLocal,particleTypes,pt_maxPerProc,&
       PRES_VAR, PRES_PART_PROP,   &
       &             POSX_PART_PROP,POSY_PART_PROP,POSZ_PART_PROP)
  call Grid_mapMeshToParticles(particles, NPART_PROPS, BLK_PART_PROP,pt_numLocal,particleTypes,pt_maxPerProc,&
       SF6_SPEC,  SF6_PART_PROP,   &
       &             POSX_PART_PROP,POSY_PART_PROP,POSZ_PART_PROP)

  !  Update derived variables by calculation
  call sim_derivedVariables(particles,pt_numLocal,particleTypes,pt_maxPerProc)

  call Timers_stop ("Particles_updateAttributes: map from mesh")

  !------------------------------------------------------------------------------

  return

end subroutine Particles_updateAttributes

