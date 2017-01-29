!!****if* source/Simulation/SimulationMain/Advect/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_init()
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get routine for initialization.
!!  Also initializes initial conditions for Advect problem.
!!
!! ARGUMENTS
!!
!!    none
!!
!! PARAMETERS
!!
!!     smallP              smallest pressure allowed
!!     smallX              smallest abundance allowed 
!!     gamma               Gamma value from the EOS
!!     sim_rhoIn           Density inside the pulse
!!     sim_rhoOut          Density outside the pulse
!!     sim_pressure        Pressure
!!     sim_velocity        Fluid velocity
!!     sim_posn            Position of the pulse center at x-axis (y=z=0)
!!     sim_width           Width of the pulse along x-axis
!!     sim_xAngle          Angle made by diaphragm normal w/x-axis (deg)
!!     sim_yAngle          Angle made by diaphragm normal w/y-axis (deg)
!!     sim_pulseFunctn     Which pulse shape function to use
!!
!!
!!***

subroutine Simulation_init()
  
  use Simulation_data
  use Driver_interface, ONLY : Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none
#include "Flash.h"
#include "constants.h"

  
  call Driver_getMype(MESH_COMM, sim_meshMe)

  
  call RuntimeParameters_get('smallp', sim_smallP)
  call RuntimeParameters_get('smallx', sim_smallX)
  
  call RuntimeParameters_get('sim_rhoin', sim_rhoIn)
  call RuntimeParameters_get('sim_rhoout', sim_rhoOut)
  call RuntimeParameters_get('gamma', sim_gamma)
  
  call RuntimeParameters_get('sim_pressure', sim_pressure)
  call RuntimeParameters_get('sim_velocity', sim_velocity)
  call RuntimeParameters_get('sim_width', sim_width)

  call RuntimeParameters_get('sim_xangle', sim_xAngle)
  call RuntimeParameters_get('sim_yangle', sim_yAngle)

  call RuntimeParameters_get('sim_posn', sim_posn)
  call RuntimeParameters_get('sim_pulseFunctn', sim_pulseFunctn)

  
  !  compute direction cosines for the pulse normal.
  sim_xAngle = sim_xAngle * 0.0174532925 ! Convert to radians.
  sim_yAngle = sim_yAngle * 0.0174532925
  
  sim_xCos = cos(sim_xAngle)
  
  if (NDIM == 1) then
     sim_xCos = 1.
     sim_yCos = 0.
     sim_zCos = 0.
     
  elseif (NDIM == 2) then
     sim_yCos = sqrt(1. - sim_xCos*sim_xCos)
     sim_zCos = 0.
     
  elseif (NDIM == 3) then
     sim_yCos = cos(sim_yAngle)
     sim_zCos = sqrt( max(0., 1. - sim_xCos*sim_xCos - sim_yCos*sim_yCos) )
  endif
end subroutine Simulation_init

