!!****if* source/physics/IncompNS/IncompNSMain/constdens/IncompNS_init
!!
!! NAME
!!
!!  IncompNS_init
!!
!!
!! SYNOPSIS
!!
!!  call IncompNS_init()
!!  
!!
!! DESCRIPTION
!! 
!!  Initialize unit scope variables which are typically the runtime parameters.
!!  This must be called once by Driver_initFlash.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!!***

subroutine IncompNS_init(restart)

  use IncompNS_data
  use ImBound_interface, ONLY : ImBound_setData
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs, &
                               Driver_getComm, Driver_getNstep, Driver_abortFlash
  use ins_interface, only : ins_getBulkVelocity
  use IncompNSstats_interface, only : IncompNSstats_init 

  use Grid_data, ONLY : gr_allPeriodic

  implicit none
  include 'Flash_mpi.h'
#include "constants.h"
#include "Flash.h"
#include "IncompNS.h"
  logical, intent(IN) :: restart

  call Driver_getMype(MESH_COMM, ins_meshMe)
  call Driver_getNumProcs(MESH_COMM, ins_meshNumProcs)
  call Driver_getComm(MESH_COMM, ins_meshComm)  

  call RuntimeParameters_get("cflflg", ins_cflflg)
  call RuntimeParameters_get("cfl", ins_cfl)
  call RuntimeParameters_get("isgs",ins_isgs)
  call RuntimeParameters_get("invRe",ins_invRe)
  call RuntimeParameters_get("sigma",ins_sigma)
  call RuntimeParameters_get("dtspec",ins_dtspec)
  call RuntimeParameters_get("intschm",ins_intschm)
  ! set scheme type flag
  select case(ins_intschm)
     case( AB2_SCHM )
        ins_intschm_type = INS_INTSCHM_MULTISTEP
     case( AB2_SCHM_V )
        ins_intschm_type = INS_INTSCHM_MULTISTEP
     case( RK3_SCHM )
        ins_intschm_type = INS_INTSCHM_RK
     case default
        call Driver_abortFlash("ins_intschm_type not known.")
  end select
  ins_vardt(:) = ins_dtspec

  call Driver_getNstep(ins_nstep)
  ins_restart=restart

  call RuntimeParameters_get("pressure_correct",ins_prescorr)
  ins_prescoeff = 0.
  if (ins_prescorr) ins_prescoeff = 1.

  call RuntimeParameters_get("vel_prolong_method",ins_prol_method)

  if (ins_meshMe .eq. MASTER_PE) then
     write(*,*) 'ins_cfl   =',ins_cfl
     write(*,*) 'ins_isgs  =',ins_isgs
     write(*,*) 'ins_invRe =',ins_invRe
     write(*,*) 'ins_sigma =',ins_sigma
     write(*,*) 'ins_dtspec=',ins_dtspec
     write(*,*) 'ins_intschm=',ins_intschm
     write(*,*) 'ins_prescoeff=',ins_prescoeff
     write(*,*) 'vel_prolong_method=',ins_prol_method
  endif

  ! Read gravity acceleration components:
  call RuntimeParameters_get("gravX",ins_gravX)
  call RuntimeParameters_get("gravY",ins_gravY)
  call RuntimeParameters_get("gravZ",ins_gravZ)
  
  ! Read pressure gradients if necessary, constant mass simulation data:
  call RuntimeParameters_get("dpdx",ins_dpdx)
  call RuntimeParameters_get("dpdy",ins_dpdy)
  call RuntimeParameters_get("dpdz",ins_dpdz)

  call RuntimeParameters_get("constantmass",ins_constmass)
  call RuntimeParameters_get("WBREF",ins_WBREF)
  call RuntimeParameters_get("area_solids",ins_area_solids)

  ! Populate old bulk velocity:
  if (ins_constmass) call ins_getBulkVelocity(ins_WBold,KAXIS)

  if(ins_meshMe==MASTER_PE) then
     if ( ins_constmass .and. (ins_WBREF .eq. 0.)) print*, 'WARNING: Constant Mass selected, but WBREF = 0.0'
     print*,'Incmp_Navier_Stokes initialized'
  endif

  ! Define value of kinematic viscosity for IB distributed forces:
  call ImBound_setData()

  ! Call initialization of ins statistics parameters
  call IncompNSstats_init()

  ! For all PERIODIC boundary conditions
  ! Shizhao Wang
  ! Nov. 07, 2014
  if (gr_allPeriodic) then
    ins_gravX = 0.
    ins_gravY = 0.
    ins_gravZ = 0.
    
!    write(*,*) 'Reset the gravity for NS solver because of the periodic boundary condition'
!    write(*,*) 'Gravity in the NS eqns.:', ins_gravX, ins_gravY, ins_gravZ
  endif

end subroutine IncompNS_init

