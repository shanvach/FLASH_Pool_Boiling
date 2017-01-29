!!****if* source/Simulation/SimulationMain/ShockCyl/sim_derivedVariables
!!
!! NAME
!!
!!  sim_derivedVariables
!!
!! SYNOPSIS
!!
!!  sim_derivedVariables
!!
!! DESCRIPTION
!!
!!   Calculates derived variables for ShockCylinder 3d problem
!!
!! ARGUMENTS
!!  
!! PARAMETERS
!!
!! NOTES
!!
!!  This routine was moved to the bottom of this file to avoid conflicts
!!    when using setup without Particles
!!
!!***

subroutine sim_derivedVariables(particles,pt_numLocal,maxParticlesPerProc)
  use Grid_interface, ONLY : Grid_getDeltas, Grid_getBlkPtr, Grid_releaseBlkPtr
  use Particles_interface, ONLY : Particles_mapFromMesh

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(IN)   :: pt_numLocal, maxParticlesPerProc
  real, intent(INOUT), dimension(NPART_PROPS,maxParticlesPerProc) :: particles
    
  integer               :: i, currentBlk

  real                  :: vx_xleft, vx_xright, vx_yleft, vx_yright, vx_zleft, vx_zright
  real                  :: vy_xleft, vy_xright, vy_yleft, vy_yright, vy_zleft, vy_zright
  real                  :: vz_xleft, vz_xright, vz_yleft, vz_yright, vz_zleft, vz_zright

  real, dimension(MDIM) :: dxArray
  real                  :: dx, dy, dz

  real                  :: twodxi, twodyi, twodzi

  real, pointer, dimension(:,:,:,:) :: solnVec

!------------------------------------------------------------------------------
!! DEV This routine shouldn't be compiled if particles aren't included
!! DEV Particles are not the default in the Config file -- same is true for Particles_updateAttributes in this directory

   do i = 1, pt_numLocal

      currentBlk=int(particles(BLK_PART_PROP,i))

      call Grid_getDeltas( currentBlk, dxArray)
      dx = dxArray(1)
      if (NDIM >= 2) then
         dy = dxArray(2)
      else
         dy = 1
      endif
      if (NDIM == 3) then
         dz = dxArray(3)
      else
         dz = 1
      endif

      twodxi = 1.e0/(2.e0*dx)
      twodyi = 1.e0/(2.e0*dy)
      twodzi = 1.e0/(2.e0*dz)

      call Grid_getBlkPtr(currentBlk,solnVec)


      call Particles_mapFromMesh(VELX_VAR, particles(POSX_PART_PROP,i) + dx,  &
                                 particles(POSY_PART_PROP,i), &
                                 particles(POSZ_PART_PROP,i), &
                                 vx_xright, currentBlk, solnVec)

      call Particles_mapFromMesh(VELX_VAR, particles(POSX_PART_PROP,i) - dx,  &
                                 particles(POSY_PART_PROP,i), &
                                 particles(POSZ_PART_PROP,i), &
                                 vx_xleft, currentBlk, solnVec)

      call Particles_mapFromMesh(VELX_VAR, particles(POSX_PART_PROP,i),  &
                                 particles(POSY_PART_PROP,i) + dy, &
                                 particles(POSZ_PART_PROP,i), &
                                 vx_yright, currentBlk, solnVec)

      call Particles_mapFromMesh(VELX_VAR, particles(POSX_PART_PROP,i),  &
                                particles(POSY_PART_PROP,i) - dy, &
                                 particles(POSZ_PART_PROP,i), &
                                 vx_yleft, currentBlk, solnVec)

      call Particles_mapFromMesh(VELX_VAR, particles(POSX_PART_PROP,i),  &
                                 particles(POSY_PART_PROP,i), &
                                 particles(POSZ_PART_PROP,i) + dz, &
                                 vx_zright, currentBlk, solnVec)

      call Particles_mapFromMesh(VELX_VAR, particles(POSX_PART_PROP,i),  &
                                 particles(POSY_PART_PROP,i), &
                                 particles(POSZ_PART_PROP,i) - dz, &
                                 vx_zleft, currentBlk, solnVec)

      call Particles_mapFromMesh(VELY_VAR, particles(POSX_PART_PROP,i) + dx,  &
                                 particles(POSY_PART_PROP,i), &
                                 particles(POSZ_PART_PROP,i), &
                                 vy_xright, currentBlk, solnVec)

      call Particles_mapFromMesh(VELY_VAR, particles(POSX_PART_PROP,i) - dx,  &
                                 particles(POSY_PART_PROP,i), &
                                 particles(POSZ_PART_PROP,i), &
                                 vy_xleft, currentBlk, solnVec)

      call Particles_mapFromMesh(VELY_VAR, particles(POSX_PART_PROP,i),  &
                                 particles(POSY_PART_PROP,i) + dy, &
                                 particles(POSZ_PART_PROP,i), &
                                 vy_yright, currentBlk, solnVec)

      call Particles_mapFromMesh(VELY_VAR, particles(POSX_PART_PROP,i),  &
                                 particles(POSY_PART_PROP,i) - dy, &
                                 particles(POSZ_PART_PROP,i), &
                                 vy_yleft, currentBlk, solnVec)

      call Particles_mapFromMesh(VELY_VAR, particles(POSX_PART_PROP,i),  &
                                 particles(POSY_PART_PROP,i), &
                                 particles(POSZ_PART_PROP,i) + dz, &
                                 vy_zright, currentBlk, solnVec)

      call Particles_mapFromMesh(VELY_VAR, particles(POSX_PART_PROP,i),  &
                                 particles(POSY_PART_PROP,i), &
                                 particles(POSZ_PART_PROP,i) - dz, &
                                 vy_zleft, currentBlk, solnVec)

      call Particles_mapFromMesh(VELZ_VAR, particles(POSX_PART_PROP,i) + dx,  &
                                 particles(POSY_PART_PROP,i), &
                                 particles(POSZ_PART_PROP,i), &
                                 vz_xright, currentBlk, solnVec)

      call Particles_mapFromMesh(VELZ_VAR, particles(POSX_PART_PROP,i) - dx,  &
                                 particles(POSY_PART_PROP,i), &
                                 particles(POSZ_PART_PROP,i), &
                                 vz_xleft, currentBlk, solnVec)

      call Particles_mapFromMesh(VELZ_VAR, particles(POSX_PART_PROP,i),  &
                                 particles(POSY_PART_PROP,i) + dy, &
                                 particles(POSZ_PART_PROP,i), &
                                 vz_yright, currentBlk, solnVec)

      call Particles_mapFromMesh(VELZ_VAR, particles(POSX_PART_PROP,i),  &
                                 particles(POSY_PART_PROP,i) - dy, &
                                 particles(POSZ_PART_PROP,i), &
                                 vz_yleft, currentBlk, solnVec)

      call Particles_mapFromMesh(VELZ_VAR, particles(POSX_PART_PROP,i),  &
                                 particles(POSY_PART_PROP,i), &
                                 particles(POSZ_PART_PROP,i) + dz, &
                                 vz_zright, currentBlk, solnVec)

      call Particles_mapFromMesh(VELZ_VAR, particles(POSX_PART_PROP,i),  &
                                 particles(POSY_PART_PROP,i), &
                                 particles(POSZ_PART_PROP,i) - dz, &
                                 vz_zleft, currentBlk, solnVec)
     call Grid_releaseBlkPtr(currentBlk,solnVec)

      particles(OMEGAX_PART_PROP,i)  = 0.5e0*((vz_yright-vz_yleft)*twodyi - (vy_zright-vy_zleft)*twodzi)
      particles(OMEGAY_PART_PROP,i)  = 0.5e0*((vx_zright-vx_zleft)*twodzi - (vz_xright-vz_xleft)*twodxi)
      particles(OMEGAZ_PART_PROP,i)  = 0.5e0*((vy_xright-vy_xleft)*twodxi - (vx_yright-vx_yleft)*twodyi)

      particles(GAMMAXY_PART_PROP,i) = ((vx_yright-vx_yleft)*twodyi + (vy_xright-vy_xleft)*twodxi)
      particles(GAMMAXZ_PART_PROP,i) = ((vx_zright-vx_zleft)*twodzi + (vz_xright-vz_xleft)*twodxi)
      particles(GAMMAYZ_PART_PROP,i) = ((vy_zright-vy_zleft)*twodzi + (vz_yright-vz_yleft)*twodyi)

      particles(DVXDX_PART_PROP,i)    = (vx_xright-vx_xleft)*twodxi
      particles(DVYDY_PART_PROP,i)    = (vy_yright-vy_yleft)*twodyi
      particles(DVZDZ_PART_PROP,i)    = (vz_zright-vz_zleft)*twodzi

   enddo

end subroutine sim_derivedVariables



