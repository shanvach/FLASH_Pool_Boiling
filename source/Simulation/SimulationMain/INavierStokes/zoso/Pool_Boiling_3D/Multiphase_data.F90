!!****if* source/physics/Multiphase/MultiphaseMain/Multiphase_data
!!
!! NAME
!!
!!  Multiphase_data
!!
!!
!! SYNOPSIS
!!
!!  MODULE Multiphase_data()
!!
!!
!! ARGUMENTS
!!
!!
!! DESCRIPTION
!!
!!  This stores data and limiter functions that are specific to the Multiphase module.
!!
!!***
 
 
module Multiphase_data

#include "Flash.h"
#include "constants.h"

  real, save :: mph_rho1
  real, save :: mph_rho2

  real, save :: mph_vis1
  real, save :: mph_vis2

  real, save :: mph_sten

  real, save :: mph_crmx, mph_crmn

  integer, save :: mph_lsit
  integer, save :: mph_inls

  integer, save :: mph_meshMe
  integer, save :: mph_meshNumProcs
  integer, save :: mph_meshComm

  real, save :: mph_thco1 ! Akash
  real, save :: mph_thco2 ! Akash
  
  real, save :: mph_cp1   ! Akash
  real, save :: mph_cp2   ! Akash

  real, save :: mph_radius

  real, save :: mph_baseRadius

  real, save :: mph_offset

  integer, save :: mph_baseCount, mph_baseCountAll

  ! For Nucleate Boiling Re-Initialization

  logical, save :: mph_isAttached 
  real, save :: mph_timeStamp

  logical, save, allocatable :: mph_isAttachedAll(:)
  logical, save, allocatable :: mph_isAttachedOld(:)
  real,    save, allocatable :: mph_timeStampAll(:)
  real,    save, allocatable :: mph_nucSiteTemp(:)

  integer, save :: mph_bcFlag(NXB,NZB,MAXBLOCKS*10)

  logical, save :: mph_redistanceFlg

  real, save    :: mph_psi(NXB+2*NGUARD,NZB+2*NGUARD*K3D,MAXBLOCKS)

  real, save    :: mph_psi_adv, mph_vlim

end module Multiphase_data
