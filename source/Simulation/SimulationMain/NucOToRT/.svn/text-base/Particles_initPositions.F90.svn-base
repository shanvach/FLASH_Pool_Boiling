!!****if* source/Simulation/SimulationMain/NucOToRT/Particles_initPositions
!!
!! NAME
!!    Particles_initPositions
!!
!! SYNOPSIS
!!
!!    Particles_initPositions( logical(inout) :: success,
!!                             logical(out)   :: updateRefine)
!!
!! DESCRIPTION
!!    Initialize particle locations.  This version is designed for the 
!!    particle based poistest problem.  The particle locations are
!!    initialized according to the density given by a 
!!    Gaussian distribution function.
!!
!! ARGUMENTS
!!
!!    success : boolean indicating whether particles positions were 
!!              successfully initialized. This is not really relevant
!!              for this version of the routine.
!! updateRefine : is true if the routine wishes to retain the already
!!                initialized particles instead of reinitializing them
!!                as the grid refines.
!!
!!
!!***

subroutine Particles_initPositions (success,updateRefine)

  use Particles_data, ONLY:  pt_numLocal, particles, pt_maxPerProc, &
       pt_attributes, &
       pt_meshVar,pt_numAttributes,&
       pt_posInitialized, pt_velInitialized, pt_logLevel, pt_meshMe
  use Simulation_data, ONLY: sim_iPartFile

#ifdef BEFORE_PARTICLEMERGE
  use Grid_interface, ONLY : Grid_moveParticlesGlobal
#else
  use Particles_interface, ONLY: Particles_updateRefinement
#endif
  use Grid_interface, ONLY : Grid_mapParticlesToMesh
  use Logfile_interface, ONLY : Logfile_stamp
  use sim_interface, ONLY : sim_readNucOutput

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Particles.h"

  logical, intent(INOUT) :: success
  logical, intent(OUT) :: updateRefine

  integer,dimension(MAXBLOCKS,1) :: particlesPerBlk

  integer       :: i, j
  integer       :: numAttrib
  integer       :: pt_numInvalid, pt_numInvBlk

  integer       :: partProp_ind
  real :: newpos(3), newVel(3)

  integer :: unused_oldLocalNumBlocks

!----------------------------------------------------------------------




  updateRefine=.false.
  if(success) return

  pt_numLocal=0

  call sim_readNucOutput(NPART_PROPS,pt_maxPerProc,pt_numLocal,iPartFile=sim_iPartFile)
  print*,'Particles_initPositions 1:local number of particles',pt_numLocal,' on',pt_meshMe
!!  print*,'Particles_initPositions 1:global number of particles',pt_numLocal,' on',pt_meshMe


  pt_numInvalid = 0
  pt_numInvBlk = 0
  do i = 1,pt_numLocal
     if (particles(TAG_PART_PROP,i) == 0) then
        pt_numInvalid = pt_numInvalid + 1
        particles(BLK_PART_PROP,i) = NONEXISTENT
     elseif (particles(BLK_PART_PROP,i) .LE. 0) then
        pt_numInvalid = pt_numInvalid + 1
        pt_numInvBlk = pt_numInvBlk + 1
     end if
     
  end do
  if (pt_numInvalid > 0) then
     print *,'[Particles_initPositions] WARNING found',pt_numInvalid,' invalid particles.'
     if (pt_logLevel > PT_LOGLEVEL_WARN_OPER) &
          call Logfile_stamp(pt_numInvalid,'[Particles_initPositions] WARNING found invalid particles')
!!$           print*,'I first have particles',pt_numLocal,NPART_PROPS,pt_meshMe
!!$           print*,'and they are',particles((/TAG_PART_PROP,BLK_PART_PROP,POSX_PART_PROP/),1:pt_numLocal)
!!$     call Grid_sortParticles(particles, NPART_PROPS, pt_numLocal, 1, pt_maxPerProc, particlesPerBlk)
!!$     pt_numLocal = pt_numLocal - pt_numInvalid
!!$     print *,pt_meshMe,': Reduced pt_numLocal by',pt_numInvalid,' to',pt_numLocal
!!$           print*,'I now have particles',pt_numLocal,NPART_PROPS,pt_meshMe
!!$           print*,'and they are',particles((/TAG_PART_PROP,BLK_PART_PROP,POSX_PART_PROP/),1:pt_numLocal)
     if (pt_numInvBlk > 0) then
        print *,' That number includes',pt_numInvBlk,' particles with valid TAG but invalid BLK attribute'
     end if
     j = 0
     do i = 1,pt_numLocal
        if ((particles(TAG_PART_PROP,i) .NE. 0) .AND. (particles(BLK_PART_PROP,i) .GT. 0)) then
           j = j + 1
           if (j < i) particles(:,j) = particles(:,i)
        end if
     end do
     pt_numLocal = j
  end if

  do i = 1,pt_numLocal
     call sim_coordTransfm(particles(POSX_PART_PROP,i),&
          particles(POSY_PART_PROP,i),&
          particles(POSZ_PART_PROP,i),&
          newpos(1), newpos(2), newpos(3),2,&
          particles(VELX_PART_PROP,i),&
          particles(VELY_PART_PROP,i),&
          particles(VELZ_PART_PROP,i),&
          newVel(1), newVel(2), newVel(3))
     particles(POSX_PART_PROP,i) = newpos(1)
     particles(POSY_PART_PROP,i) = newpos(2)
     particles(POSZ_PART_PROP,i) = newpos(3)
     particles(VELX_PART_PROP,i) = newVel(1)
     particles(VELY_PART_PROP,i) = newVel(2)
     particles(VELZ_PART_PROP,i) = newVel(3)
  end do


  success = .true.
  pt_posInitialized = success

#ifdef BEFORE_PARTICLEMERGE
  call Grid_moveParticlesGlobal(particles,pt_numLocal,pt_maxPerProc)
#else
  call Particles_updateRefinement(unused_oldLocalNumBlocks)
#endif
  print*,'Particles_initPositions 2:local number of particles',pt_numLocal,' on',pt_meshMe

  numAttrib = 0
  do i = 1,pt_numAttributes
     if(pt_meshVar(PT_MAP,i)==PARTICLEMAP_UNK) then
        numAttrib=numAttrib+1
        partProp_ind = pt_attributes(i)
        ! Prevent initialization of particle velocities from the mesh in Particles_initData,
        ! since we want to go the other way and later initialize mesh velocities from the
        ! particle velocity attributes. - KW
        if (partProp_ind==VELX_PART_PROP) pt_velInitialized = .TRUE.
     end if
  end do

  print*,'Particles_initPositions 3:local number of particles',pt_numLocal,' on',pt_meshMe
  if (pt_meshMe.EQ.MASTER_PE) print*,'the local number of particles is ',pt_numLocal


  return

end subroutine Particles_initPositions
