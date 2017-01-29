!!****if* source/Simulation/SimulationMain/VParticles_DPD/Simulation_init
!!
!! NAME
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!  Simulation_init( )
!!
!! DESCRIPTION   
!!   Initialize all the runtime parameters needed for testing virtual 
!!   particles functionality
!!
!! ARGUMENTS
!!
!!   
!!
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Particles_data, ONLY:  pt_numLocal, particles, pt_maxPerProc, &
       pt_xmin, pt_ymin, pt_zmin,pt_xmax, pt_ymax, pt_zmax,&
       pt_posAttrib,pt_velNumAttrib, pt_velAttrib,pt_typeInfo, pt_meshMe
  implicit none

#include "constants.h"
#include "Flash.h"
  
  
  integer :: i, j, status
 

  sig = 3.
  rc  = 1.;
  lo  = 0.25;
  ks  = 10.;
  kd  = 5.;
  !k_bolz=1.38e-21;
  kt=1.
  sig_sq= sig*sig;   
  lambda=0.5;
  call machineepsilon 
  call RuntimeParameters_get("pt_NumPart",   pt_NumPart)
  write(*,*)'The number of particles defined in par file= ',pt_NumPart
  call RuntimeParameters_get("pt_numBodies", pt_NumBodies)
 write(*,*)'The number of bodies defined in par file= ',pt_NumBodies
  call RuntimeParameters_get("pt_numBTypes",pt_NumBTypes)
 write(*,*)'The number of beadtypes defined in par file= ',pt_NumBTypes
!!$  call RuntimeParameters_get("sim_vz_pert",   sim_vz_pert )
 call RuntimeParameters_get("pt_BodyTypes",pt_BodyTypes)
 write(*,*)'The number of bodytypes defined in par file= ',pt_BodyTypes
 !call RuntimeParameters_get("sigma",sig)         ! noise amplitude 
end subroutine Simulation_init
