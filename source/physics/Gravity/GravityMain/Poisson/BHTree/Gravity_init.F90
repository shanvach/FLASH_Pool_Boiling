!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_init
!!
!! NAME
!!
!!  Gravity_init
!!
!! 
!! SYNOPSIS
!!
!!  call Gravity_init()
!!
!!
!! DESCRIPTION
!!
!!  Initialize the Barnes-Hut tree Poisson solver.  Read in any of the
!!  runtime parameters for this solver.  All solver common data
!!  is stored in the tree_common module
!!
!!***

subroutine Gravity_init()

  use Gravity_data, ONLY : newton, grav_boundary, grav_geometry,&
       grav_poisfact, grv_meshMe, grv_meshNumProcs, grv_meshComm,&
       updateGravity, useGravity
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype,&
       Driver_getComm, Driver_getNumProcs
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
    RuntimeParameters_mapStrToInt
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  implicit none
#include "constants.h"

  character(len=MAX_STRING_LENGTH) :: strGeometry
  character(len=MAX_STRING_LENGTH) :: grav_boundary_type, grav_boundary_type_x,&
     grav_boundary_type_y, grav_boundary_type_z

  call Driver_getMype(MESH_COMM,grv_meshMe)
  call Driver_getComm(MESH_COMM,grv_meshComm)
  call Driver_getNumProcs(MESH_COMM,grv_meshNumProcs)

  call RuntimeParameters_get("geometry", strGeometry)
  call RuntimeParameters_mapStrToInt(strGeometry, grav_geometry)

  

  call RuntimeParameters_get("useGravity", useGravity)
  call RuntimeParameters_get("updateGravity", updateGravity)

  call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)
  call RuntimeParameters_get("grav_boundary_type_x", grav_boundary_type_x)
  call RuntimeParameters_get("grav_boundary_type_y", grav_boundary_type_y)
  call RuntimeParameters_get("grav_boundary_type_z", grav_boundary_type_z)
  
  call PhysicalConstants_get("Newton", newton)


  !! FUTURE : need to add scale factor for cosmology
  !grav_scaleFactor = 1
  grav_poisfact = 4. * PI * Newton

  select case (grav_boundary_type)
    case("isolated")
      grav_boundary = ISOLATED
    case("periodic")
      grav_boundary = PERIODIC
    case("mixed")
      select case (grav_boundary_type_x)
        case("isolated")
          grav_boundary(1) = ISOLATED
          grav_boundary(2) = ISOLATED
        case("periodic")
          grav_boundary(1) = PERIODIC
          grav_boundary(2) = PERIODIC
        case default
          call Driver_abortFlash("Gravity_init: unrecognized or unsupported gravity boundary type in x direction")
      end select
      select case (grav_boundary_type_y)
        case("isolated")
          grav_boundary(3) = ISOLATED
          grav_boundary(4) = ISOLATED
        case("periodic")
          grav_boundary(3) = PERIODIC
          grav_boundary(4) = PERIODIC
        case default
          call Driver_abortFlash("Gravity_init: unrecognized or unsupported gravity boundary type in y direction")
      end select
      select case (grav_boundary_type_z)
        case("isolated")
          grav_boundary(5) = ISOLATED
          grav_boundary(6) = ISOLATED
        case("periodic")
          grav_boundary(5) = PERIODIC
          grav_boundary(6) = PERIODIC
        case default
          call Driver_abortFlash("Gravity_init: unrecognized or unsupported gravity boundary type in z direction")
      end select
    case default
      call Driver_abortFlash("Gravity_init: unrecognized or unsupported gravity boundary type")
  end select

  return
end subroutine Gravity_init
