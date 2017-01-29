!!****if* source/physics/sourceTerms/Heat/HeatMain/StatPlusGauss/Heat_init
!!
!! NAME
!!  
!!  Heat_init
!!
!!
!! SYNOPSIS
!! 
!!  call Heat_init()
!!
!!  
!! DESCRIPTION
!!
!!  Perform various initializations (apart from the problem-dependent ones)
!!  for the heat module.
!!
!!
!! ARGUMENTS
!!
!!   
!!
!!***
subroutine Heat_init()

  use Heat_data, ONLY: ht_x0, ht_y0, ht_z0,&
       ht_stat, ht_q, ht_sig, ht_tstar, ht_t0, ht_tau, ht_tmin, &
       useHeat, ht_meshMe, ht_numProcs

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_getMype,Driver_getNumProcs

  implicit none

#include "constants.h"

  ! Everybody should know theseb

  call Driver_getMype(MESH_COMM,ht_meshMe)
  call Driver_getNumProcs(MESH_COMM, ht_numProcs)


  call RuntimeParameters_get("statheat", ht_stat)
  call RuntimeParameters_get("qheat", ht_q)
  call RuntimeParameters_get("x0heat", ht_x0)
  call RuntimeParameters_get("y0heat", ht_y0)
  call RuntimeParameters_get("z0heat", ht_z0)
  call RuntimeParameters_get("sigheat", ht_sig)
  call RuntimeParameters_get("tstar", ht_tstar)
  call RuntimeParameters_get("t0heat", ht_t0)
  call RuntimeParameters_get( "tau", ht_tau)
  call RuntimeParameters_get("theatmin", ht_tmin)
  call RuntimeParameters_get("useHeat",useHeat)
  if (.not. useHeat) then
     write(6,*)'WARNING:  You have included the Heat unit but have set '
     write(6,*)'   the runtime parameter useHeat to FALSE'
     write(6,*)'   No heating will occur but Heat_init will continue.'
  end if




  return
end subroutine Heat_init
