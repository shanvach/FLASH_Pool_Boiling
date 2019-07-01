!!****if* source/physics/Multiphase/MultiphaseMain/Multiphase_init
!!
!! NAME
!!
!!  Multiphase_init
!!
!!
!! SYNOPSIS
!!
!!  call Multiphase_init()
!!  
!!
!! DESCRIPTION
!! 
!!
!!***

subroutine Multiphase_init( )

  use Multiphase_data, ONLY : mph_rho1,mph_rho2,mph_sten, &
                              mph_vis1,mph_vis2,mph_lsit, mph_inls, &
                              mph_meshMe, mph_meshNumProcs, mph_meshComm, &
                              mph_jet_tstamp, mph_jet_twait, mph_jet_period, &
                              mph_dist_flag, mph_vel_scalar
 
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs, &
                               Driver_getComm

  use Driver_data, only: dr_simTime

  implicit none
  include 'Flash_mpi.h'
#include "constants.h"
#include "Flash.h"

  call Driver_getMype(MESH_COMM, mph_meshMe)
  call Driver_getNumProcs(MESH_COMM, mph_meshNumProcs)
  call Driver_getComm(MESH_COMM, mph_meshComm)


  call RuntimeParameters_get("rho1",mph_rho1)
  call RuntimeParameters_get("rho2",mph_rho2)
  call RuntimeParameters_get("vis1",mph_vis1)
  call RuntimeParameters_get("vis2",mph_vis2)
  call RuntimeParameters_get("sten",mph_sten)
  call RuntimeParameters_get("lsit",mph_lsit)
  call RuntimeParameters_get("inls",mph_inls)

 if (mph_meshMe .eq. MASTER_PE) then
     write(*,*) 'mph_rho1=',mph_rho1
     write(*,*) 'mph_rho2=',mph_rho2
     write(*,*) 'mph_vis1=',mph_vis1
     write(*,*) 'mph_vis2=',mph_vis2
     write(*,*) 'mph_sten=',mph_sten
     write(*,*) 'mph_lsit=',mph_lsit
     write(*,*) 'mph_inls=',mph_inls
  endif

  mph_dist_flag  = .true.
  mph_jet_tstamp = dr_simTime
  mph_jet_twait  = 50.0
  mph_jet_period = 50.0
  mph_vel_scalar = 1.0

end subroutine Multiphase_init
