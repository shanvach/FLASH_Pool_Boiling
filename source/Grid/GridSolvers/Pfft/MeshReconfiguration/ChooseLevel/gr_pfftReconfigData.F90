!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/ChooseLevel/gr_pfftReconfigData
!!
!! NAME
!!  gr_pfftReconfigData
!!
!! SYNOPSIS
!!  use gr_pfftReconfigData
!!
!!  This includes subunit scope data for the pfft data movement
!!  
!!***

module gr_pfftReconfigData
#include "constants.h"
#include "Flash.h"
  use gr_pfftMessageList, ONLY : message_list
  implicit none

  !We have devised pfft_procLookup which is an array of pointers to a 
  !a small datatype consisting of 2 integers. It allows us to extract 
  !information about a PFFT processor and its neighbors directly 
  !(i.e. no extensive search required, as is the case with pure arrays).

  !Example usage:
  !pfft_procLookup(JAXIS) % procInfo(3) % globalStartGridPoint
  !corresponds to information for the 4th (it is zero based) 
  !processor along the JAXIS.
  type TProcInfo
     integer :: globalStartGridPoint
     integer :: globalEndGridPoint
  end type TProcInfo

  type TPtrProcInfo
     type(TProcInfo), dimension(:), pointer :: procInfo
  end type TPtrProcInfo

  type(TPtrProcInfo), save, dimension(1:NDIM) :: pfft_procLookup


  integer, save :: pfft_pencilSize !Size of pencil assigned to pfft_myPE.

  logical, save :: pfft_debugMode   !Debug mode: Add debugging checks.
  logical, save :: pfft_logMode !Log mode: Add process specific logging.
  integer, save :: pfft_logUnit !Log file unit number.

  !Integers representing FLASH grid unk variables.
  integer, save :: pfft_srcGridVar, pfft_solnGridVar

  !Pointer to a particular 1D PFFT array:
  real, save, dimension(:), pointer :: pfft_pfftBuf

  !There are 2 main lists in the code:
  type(message_list), save, pointer :: pfft_flashMsgList, pfft_pencilMsgList

  integer, save :: pfft_mpiMetadataType
  integer, save :: pfft_minRefLev, pfft_maxRefLev, pfft_oneRefLev
  integer, parameter :: pfft_mpiTag1 = 123
  integer, parameter :: pfft_mpiTag2 = 124
  integer, parameter :: pfft_mpiTag3 = 125
  integer, parameter :: pfft_mpiTag4 = 126
  integer, parameter :: pfft_mpiTag5 = 127

end module gr_pfftReconfigData
