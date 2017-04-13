!!  Multiphase_interface
!!
!! SYNOPSIS
!!
!!  use Multiphase_interface
!!
!! DESCRIPTION
!!
!! This is the header file for the Multiphase module
!! module that defines its public interfaces.
!!
!!***

Module Multiphase_interface

  implicit none

#include "constants.h"
#include "Flash.h"


  interface
  subroutine Multiphase_init()
  implicit none
  end subroutine
  end interface

  interface

  subroutine Multiphase(blockCount,blockList,timeEndAdv,dt,dtOld,sweepOrder,mph_flag)
  implicit none
  integer, INTENT(IN) :: sweepOrder
  integer, INTENT(INOUT) :: blockCount
  integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList !blockCount
  real,    INTENT(IN) :: timeEndAdv,dt,dtOld
  integer, intent(in) :: mph_flag
  end subroutine
  end interface

end module Multiphase_interface
