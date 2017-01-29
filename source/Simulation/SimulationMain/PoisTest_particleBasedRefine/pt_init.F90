!!****if* source/Simulation/SimulationMain/PoisTest_particleBasedRefine/pt_init
!!
!! NAME
!!    pt_init
!!
!! SYNOPSIS
!!
!!    pt_init()
!!
!! DESCRIPTION
!!
!!    This routine puts values in the runtime parameters and other 
!!    ParticlesInitialization subunit scope variables for the 
!!    particle based poistest problem.
!!
!!
!! PARAMETERS
!!
!!    pt_numX:      number of particles along physical x-axis of domain
!!    pt_numY:      number of particles along physical y-axis of domain
!!    pt_numZ:      number of particles along physical z-axis of domain
!!
!!
!! NOTES
!!
!!***

subroutine pt_init()


  use pt_initData, ONLY:  pt_numX, pt_numY, pt_numZ
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  
  implicit none
 
  ! Initial distribution parameters.
  call RuntimeParameters_get ("pt_numX", pt_numX)
  call RuntimeParameters_get ("pt_numY", pt_numY)
  call RuntimeParameters_get ("pt_numZ", pt_numZ)
  
  return
  
end subroutine pt_init
