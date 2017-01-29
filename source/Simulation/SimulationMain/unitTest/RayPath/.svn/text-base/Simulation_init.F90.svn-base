!!****if* source/Simulation/SimulationMain/unitTest/RayPath/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data, for the 
!!  unit test tracing the path of all rays incident on the domain
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Initializes initial incidence point and angle for all rays 
!!  traversing through the domain
!!
!! ARGUMENTS
!!
!!   
!!
!! PARAMETERS
!!
!!  sim_ilBnd         The lower IAXIS bound of the refracting slab
!!  sim_iubnd         The upper IAXIS bound 
!!  sim_jlBnd         The lower JAXIS bound of the refracting slab
!!  sim_jubnd         The upper JAXIS bound 
!!  sim_klBnd         The lower KAXIS bound of the refracting slab
!!  sim_kubnd         The upper KAXIS bound 
!!  sim_refract       Refractive index of the slab
!!  sim_refractType   the type of refraction (constant, or with constant
!!                    gradient)
!!  sim_numRay        The number of rays 
!!  sim_rayIncidence  array containing the initial entry point and incidence
!!                    angles of each ray.
!!  sim_fileRay       filename of the file containing the incidence info
!!
!!***

subroutine Simulation_init()

  use Simulation_data , ONLY : sim_ilBnd, sim_iuBnd, sim_jlBnd, sim_juBnd,&
                               sim_klBnd, sim_kuBnd, sim_refractType,&
                               sim_refract, sim_numRay, sim_rayIncidence,&
                               sim_fileRay, sim_slabBndBox
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "constants.h"
#include "Flash.h"

  

  integer :: i

  call RuntimeParameters_get('sim_ilBnd', sim_ilBnd)
  call RuntimeParameters_get('sim_iuBnd', sim_iuBnd)
  call RuntimeParameters_get('sim_jlBnd', sim_jlBnd)
  call RuntimeParameters_get('sim_juBnd', sim_juBnd)
  call RuntimeParameters_get('sim_refract',sim_refract)
  call RuntimeParameters_get('sim_refractType',sim_refractType)
  call RuntimeParameters_get('sim_rayNum', sim_rayNum)
  call RuntimeParameters_get('sim_fileRay',sim_fileRay)
  call RuntimeParameters_get('sim_domainBndBox(LOW,IAXIS),'xmin')
  call RuntimeParameters_get('sim_domainBndBox(HIGH,IAXIS),'xmax')
  call RuntimeParameters_get('sim_domainBndBox(LOW,JAXIS),'ymin')
  call RuntimeParameters_get('sim_domainBndBox(HIGH,JAXIS),'ymax')
  call RuntimeParameters_get('sim_domainBndBox(LOW,KAXIS),'zmin')
  call RuntimeParameters_get('sim_domainBndBox(HIGH,KAXIS),'zmax')

  allocate(sim_rayIncidence(IAXIS:KAXIS+2,sim_rayNum)
  allocate(sim_rayFinal(IAXIS:KAXIS+1,sim_rayNum)
  open(1,file=sim_fileRay)
  do i = 1,sim_rayNum
     !! Read in the initial position of incidence, and the slope
     !! of the ray along x-y and x-z axes
     read(1,*)sim_rayIncidence(IAXIS:KAXIS+2,i)
  end do
  sim_slabBndBox(LOW,IAXIS)=sim_ilBnd
  sim_slabBndBox(HIGH,IAXIS)=sim_iuBnd
  sim_slabBndBox(LOW,JAXIS)=sim_jlBnd
  sim_slabBndBox(HIGH,JAXIS)=sim_juBnd
  sim_slabBndBox(LOW,KAXIS)=sim_klBnd
  sim_slabBndBox(HIGH,KAXIS)=sim_kuBnd

end subroutine Simulation_init
