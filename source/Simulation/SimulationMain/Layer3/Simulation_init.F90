!!****if* source/Simulation/SimulationMain/Layer3/Simulation_init
!!
!! NAME
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!  call Simulation_init( )
!!
!!
!! DESCRIPTION
!!  Initializes all the parameters needed for the 3-Layer Problem
!!  of Remington et al.
!!
!! ARGUMENTS
!!
!!  none
!!
!!***

subroutine Simulation_init()
  use Simulation_data, ONLY : sim_xmin,sim_xmax,sim_ymin,sim_ymax,&
                              sim_small,sim_pi,sim_xzn1d,sim_1dModel,&
                              sim_gamcu,sim_gamcf,sim_gamch, sim_meshMe
  use Driver_interface, ONLY : Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Multispecies_interface, ONLY : Multispecies_getProperty
  implicit none
#include "constants.h"
#include "Flash.h"

  
  call Driver_getMype(MESH_COMM, sim_meshMe)


  !! These parameters are defined by other units
  call RuntimeParameters_get( "xmin", sim_xmin)
  call RuntimeParameters_get( "xmax", sim_xmax)
  call RuntimeParameters_get( "ymin", sim_ymin)
  call RuntimeParameters_get( "ymax", sim_ymax)
  call RuntimeParameters_get( "small",sim_small)

  !! Get the value of pi
  sim_pi = PI

  !! Get the gamma values for the three materials
  call Multispecies_getProperty(CU_SPEC,GAMMA,sim_gamcu)
  call Multispecies_getProperty(CH_SPEC,GAMMA,sim_gamch)
  call Multispecies_getProperty(CF_SPEC,GAMMA,sim_gamcf)

  !! initialize the 1-d model
  print *,'calling init_1d'

  call sim_init1d(N1D,NUNK_VARS,sim_xzn1d,sim_1dModel)
  if (sim_meshMe == MASTER_PE) then
     open(unit=79,file='test2a.out',status='unknown')
     write(6,*)'n1d = ',N1D
     do i = 1,N1D
        write(79,102)sim_xzn1d(i),sim_1dModel(i,1:NPROP_VARS)
     enddo
     close(79)
  endif
  write(6,*)'called init_1d',sim_meshMe

end subroutine Simulation_init
