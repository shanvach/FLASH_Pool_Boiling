!!****if* source/Simulation/SimulationMain/Advect/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!  use Simulation_data
!!
!! DESCRIPTION
!!
!!  Stores the local data for Simulation setup: Advect
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
!!***

module Simulation_data

  implicit none
  
  integer, save :: sim_meshMe

  !! *** Runtime Parameters *** !!
  
  real, save :: sim_rhoIn, sim_rhoOut, sim_pressure, sim_velocity
  real, save :: sim_width, sim_xAngle, sim_yAngle, sim_posn
  real, save :: sim_gamma, sim_smallP, sim_smallX
  integer, save :: sim_pulseFunctn

  !! *** Variables pertaining to this Simulation *** !!
  real, save :: sim_xCos, sim_yCos, sim_zCos

end module Simulation_data
