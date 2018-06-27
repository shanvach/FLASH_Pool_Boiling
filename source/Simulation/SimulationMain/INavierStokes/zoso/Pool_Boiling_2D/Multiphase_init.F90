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

subroutine Multiphase_init()

  use Multiphase_data, ONLY : mph_rho1,mph_rho2,mph_sten, &
                              mph_vis1,mph_vis2,mph_lsit, mph_inls, &
                              mph_meshMe, mph_meshNumProcs, mph_meshComm, &
                              mph_thco1,mph_thco2,mph_cp1,mph_cp2,mph_isAttached,mph_isAttachedAll,mph_isAttachedOld, &
                              mph_timeStampAll,mph_nucSiteTemp,mph_vlim,mph_psi_adv ! Akash
 

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs, &
                               Driver_getComm

  use Driver_data, ONLY: dr_restart

  use Simulation_data, ONLY: sim_nucSiteDens

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
  call RuntimeParameters_get("thco1",mph_thco1)
  call RuntimeParameters_get("thco2",mph_thco2)
  call RuntimeParameters_get("cp1",mph_cp1)
  call RuntimeParameters_get("cp2",mph_cp2)

  !mph_cp1 = mph_cp1*mph_rho1
  !mph_cp2 = mph_cp2*mph_rho2

 if (mph_meshMe .eq. MASTER_PE) then
     write(*,*) 'mph_rho1=',mph_rho1
     write(*,*) 'mph_rho2=',mph_rho2
     write(*,*) 'mph_vis1=',mph_vis1
     write(*,*) 'mph_vis2=',mph_vis2
     write(*,*) 'mph_sten=',mph_sten
     write(*,*) 'mph_lsit=',mph_lsit
     write(*,*) 'mph_inls=',mph_inls
     write(*,*) 'mph_thco1=',mph_thco1
     write(*,*) 'mph_thco2=',mph_thco2
     write(*,*) 'mph_cp1=',mph_cp1
     write(*,*) 'mph_cp2=',mph_cp2
     write(*,*) 'sim_nucSiteDens=',sim_nucSiteDens
  endif

  mph_vlim     = 0.2
  mph_psi_adv  = (90.0/180.0)*acos(-1.0)

  if(dr_restart .eqv. .FALSE.) then

  allocate(mph_timeStampAll(sim_nucSiteDens))
  allocate(mph_isAttachedAll(sim_nucSiteDens))
  allocate(mph_isAttachedOld(sim_nucSiteDens))
  allocate(mph_nucSiteTemp(sim_nucSiteDens))

  mph_isAttached       = .false.
  mph_isAttachedAll(:) = .false.
  mph_isAttachedOld(:) = .false.
  mph_timeStampAll(:)  = 0.0
  mph_nucSiteTemp(:)   = 0.0

  if(sim_nucSiteDens .gt. 20)  mph_timeStampAll(21:30)   = 0.2
  if(sim_nucSiteDens .gt. 30)  mph_timeStampAll(31:40)   = 0.4
  if(sim_nucSiteDens .gt. 40)  mph_timeStampAll(41:50)   = 0.6
  if(sim_nucSiteDens .gt. 50)  mph_timeStampAll(51:60)   = 0.8
  if(sim_nucSiteDens .gt. 60)  mph_timeStampAll(61:70)   = 1.0
  if(sim_nucSiteDens .gt. 70)  mph_timeStampAll(71:80)   = 1.2
  if(sim_nucSiteDens .gt. 80)  mph_timeStampAll(81:90)   = 1.4
  if(sim_nucSiteDens .gt. 90)  mph_timeStampAll(91:100)  = 1.6
  if(sim_nucSiteDens .gt. 100) mph_timeStampAll(101:110) = 1.8
  if(sim_nucSiteDens .gt. 110) mph_timeStampAll(111:120) = 2.0
  if(sim_nucSiteDens .gt. 120) mph_timeStampAll(121:130) = 2.2
  if(sim_nucSiteDens .gt. 130) mph_timeStampAll(131:139) = 2.4

  else

  allocate(mph_isAttachedAll(sim_nucSiteDens))
  allocate(mph_isAttachedOld(sim_nucSiteDens))
  allocate(mph_nucSiteTemp(sim_nucSiteDens))

  mph_isAttached       = .false.
  mph_isAttachedAll(:) = .false.
  mph_isAttachedOld(:) = .false.
  mph_nucSiteTemp(:)   = 0.0

  end if

end subroutine Multiphase_init
