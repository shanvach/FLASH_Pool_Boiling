!!****if* source/Simulation/SimulationMain/AdvectMassScalars/Simulation_init
!!
!! NAME
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!  call Simulation_init()
!!
!! DESCRIPTION   
!!     Initialize all the parameters needed for the advection of mass scalars simulation
!!
!! ARGUMENTS
!!      None.  All data passed through Simulation_data
!!
!! PARAMETERS   
!!      Described in the Config file
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use Driver_interface, ONLY : Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
  
  
  integer :: i, j, status
  real :: deg2rad

  call Driver_getMype(MESH_COMM, sim_meshMe)

  call RuntimeParameters_get("xmin", sim_xmin)
  call RuntimeParameters_get("xmax", sim_xmax)

  call RuntimeParameters_get("ymin", sim_ymin)
  call RuntimeParameters_get("ymax", sim_ymax)

  call RuntimeParameters_get("zmin", sim_zmin)
  call RuntimeParameters_get("zmax", sim_zmax)

  call RuntimeParameters_get("planar", sim_planar)

  call RuntimeParameters_get("smallp", sim_smallp)
  call RuntimeParameters_get("smallp", sim_smallp)
  call RuntimeParameters_get("smallp", sim_smallp)
  call RuntimeParameters_get("smallp", sim_smallp)
  call RuntimeParameters_get("smallp", sim_smallp)
  call RuntimeParameters_get("smallp", sim_smallp)

  call RuntimeParameters_get("smallp", sim_smallp)
  call RuntimeParameters_get("smallx", sim_smallx)
  call RuntimeParameters_get("gamma",  sim_gamma)
  call RuntimeParameters_get("rhoin",  sim_rhoin)
  call RuntimeParameters_get("rhoout", sim_rhoout)
  call RuntimeParameters_get("msin",   sim_msin)
  call RuntimeParameters_get("msout",  sim_msout)
  call RuntimeParameters_get("pressure", sim_pressure)
  call RuntimeParameters_get("velocity", sim_velocity)
  call RuntimeParameters_get("width",  sim_width)
  call RuntimeParameters_get("phase",  sim_phase)
  call RuntimeParameters_get("xangle", sim_xangle)
  call RuntimeParameters_get("yangle", sim_yangle)
  call RuntimeParameters_get("zangle", sim_zangle)
  call RuntimeParameters_get("posn",   sim_posn)

  call RuntimeParameters_get("pulse_fctn", sim_pulse_fctn)
  call RuntimeParameters_get("pulse_fctn_ms1", sim_pulse_fctn_ms1)
  call RuntimeParameters_get("pulse_fctn_ms2", sim_pulse_fctn_ms2)
  call RuntimeParameters_get("pulse_fctn_ms3", sim_pulse_fctn_ms3)
  call RuntimeParameters_get("pulse_fctn_ms4", sim_pulse_fctn_ms4)
  call RuntimeParameters_get("pulse_fctn_ms5", sim_pulse_fctn_ms5)

  sim_pi = PI

  deg2rad = sim_pi/180.e0

  sim_xangle = sim_xangle * deg2rad ! Convert to radians.
  sim_yangle = sim_yangle * deg2rad
  sim_zangle = sim_zangle * deg2rad

  sim_xcos = cos(sim_xangle)

  if (NDIM == 1) then
     sim_xcos = 1.e0
     sim_ycos = 0.e0
     sim_zcos = 0.e0

  elseif (NDIM == 2) then
     sim_ycos = sqrt(1.e0 - sim_xcos**2)
     sim_zcos = 0.e0

  elseif (NDIM == 3) then
     sim_ycos = cos(sim_yangle)
     sim_zcos = cos(sim_zangle) !sqrt( max(0.e0, 1.e0 - sim_xcos**2 - sim_ycos**2) )
  endif

  if ( NDIM < 2 ) sim_planar = .true.
  
  if ( sim_planar ) then
     sim_multid = 0.e0
  else
     sim_multid = 1.e0
  end if

  sim_twoPi = 2.e0 * sim_pi
  
  sim_xcent = sim_xmin + sim_posn*(sim_xmax - sim_xmin)
  sim_ycent = sim_ymin + sim_posn*(sim_ymax - sim_ymin)
  sim_zcent = sim_zmin + sim_posn*(sim_zmax - sim_zmin)

  ! write out to stdout listing the parameters
  if (sim_meshMe == MASTER_PE) then

     write (*,*)
     call Logfile_stamp( & 
          &           "initializing for planar advection problem.", '[Simulation_init]')
     write (*,*)  & 
          &           'flash:  initializing for mass scalar planar advection problem.'
     write (*,*)
     
     write (*,1)  & 
          &           'xcos  = ', sim_xcos,  & 
          &           'ycos  = ', sim_ycos,  & 
          &           'zcos  = ', sim_zcos

     write (*,1) 'rhoin = ', sim_rhoin,  & 
             &                  'rhoout= ', sim_rhoout,  & 
             &                  'press = ', sim_pressure

     write (*,1) 'veloc = ', sim_velocity,  & 
          &                  'width = ', sim_width, &
          &                  'phase = ', sim_phase

     write (*,2) 'posn  = ', sim_posn,  & 
          &                  'gamma = ', sim_gamma,  & 
          &                  'ndim  = ', NDIM

     write (*,*)

1    format (1X, 4(a8, es13.7, :, 1X))
2    format (1X, 2(a8, es13.7, 1X), A8, I13)

  end if

#ifdef MS1_MSCALAR
  sim_ims1 = MS1_MSCALAR
#else
  sim_ims1 = 0
#endif


#ifdef MS2_MSCALAR
  sim_ims2 = MS2_MSCALAR
#else
  sim_ims2 = 0
#endif


#ifdef MS3_MSCALAR
  sim_ims3 = MS3_MSCALAR
#else
  sim_ims3 = 0
#endif


#ifdef MS4_MSCALAR
  sim_ims4 = MS4_MSCALAR
#else
  sim_ims4 = 0
#endif


#ifdef MS5_MSCALAR
  sim_ims5 = MS5_MSCALAR
#else
  sim_ims5 = 0
#endif


  return
end subroutine Simulation_init
