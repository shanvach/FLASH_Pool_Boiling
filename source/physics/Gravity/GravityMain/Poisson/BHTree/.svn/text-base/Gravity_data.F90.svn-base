!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_data
!!
!! NAME
!!
!!  Gravity_data
!!  
!! SYNOPSIS
!!
!!  use Gravity_data
!!
!! DESCRIPTION
!!
!!  This modules stores the data for the Gravity unit
!!
!!***

module Gravity_data

#include "constants.h"

  character(len=MAX_STRING_LENGTH), save :: grav_boundary_type !string boundary condition
  integer, save :: grav_boundary(6)  !integer boundary condition

  integer, save :: grav_geometry  !mesh geometry
  integer, save :: grv_meshMe, grv_meshNumProcs, grv_meshComm
  integer, save :: grv_commSize=1

  logical, save :: useGravity, updateGravity
  logical, save :: grav_temporal_extrp !extrapolate or otherwise rescale
  logical, save :: grav_unjunkPden

  real,    save :: grav_poisfact
  real,    save :: newton

end module Gravity_data
