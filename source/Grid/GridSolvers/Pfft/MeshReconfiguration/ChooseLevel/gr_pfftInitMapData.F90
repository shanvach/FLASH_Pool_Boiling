!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/ChooseLevel/gr_pfftInitMapData
!!
!! NAME 
!!
!! gr_pfftInitMapData
!!
!! SYNOPSIS
!!
!! gr_pfftInitMapData()
!!
!! DESCRIPTION 
!!
!! Intialises variables, arrays and data structure which will be used
!! when we map from FLASH grid -> Pencil grid and Pencil grid -> FLASH grid. 
!!
!! ARGUMENTS
!!
!! SIDE EFFECTS 
!!
!! Initialisation of module level variables, arrays, lists:
!!
!! NOTES
!! 
!! We construct an MPI datatype in this routine which describes the metadata 
!! that we exchange between nodes.
!!
!!***

#define LOG_PFFT

subroutine gr_pfftInitMapData()
#include "constants.h"
#include "Flash.h"
  use Logfile_interface, ONLY : Logfile_open
  use gr_pfftReconfigData, ONLY : pfft_srcGridVar, pfft_solnGridVar, &
       pfft_flashMsgList, pfft_pencilMsgList, pfft_debugMode, pfft_logMode, &
       pfft_oneRefLev, pfft_minRefLev, pfft_maxRefLev, &
       pfft_mpiMetadataType, pfft_logUnit
  use gr_pfftData, ONLY : pfft_myPE, pfft_numProcs, pfft_comm
  use gr_pfftMessageList, ONLY : initialise_list
  use Grid_data, ONLY : gr_oneRefLev
#if defined(FLASH_GRID_PARAMESH)
  use tree, ONLY : lrefine_min, lrefine_max
#endif
  use gr_pfftMessageNode, ONLY : &
       m_logMessageNodeEvents => m_logEvents, &
       m_logMessageUnitNum => m_unitNum, &
       m_countMessageNodes => m_countID
  use gr_pfftFragmentNode, ONLY : &
       m_logFragmentNodeEvents => m_logEvents, &
       m_logFragmentUnitNum => m_unitNum, &
       m_countFragmentNodes => m_countID
  use Driver_interface, ONLY : Driver_checkMPIErrorCode
  implicit none
  include "Flash_mpi.h"
  integer, dimension(1) :: oldTypes, blockCounts, offsets
  integer :: ierr

  !Only allow pre-processor definitions in this file and gr_pfftFreeMap.
  !This minimises the amount of dead code which is a problem in other 
  !parts of FLASH.
#ifdef DEBUG_PFFT
  pfft_debugMode = .true.
#else
  pfft_debugMode = .false.
#endif

#ifdef LOG_PFFT
  pfft_logMode = .true.
  call Logfile_open(pfft_logUnit,.true.)
#else
  pfft_debugMode = .false.
#endif

#if defined(FLASH_GRID_UG)
  pfft_minRefLev = 1; pfft_maxRefLev = 1
#elif defined(FLASH_GRID_PARAMESH)
  pfft_minRefLev = lrefine_min; pfft_maxRefLev = lrefine_max
#endif
  pfft_srcGridVar = NONEXISTENT; pfft_solnGridVar = NONEXISTENT
  pfft_oneRefLev = gr_oneRefLev


  !Initialise module level variables taking charge of logging events.
  m_logMessageNodeEvents = pfft_logMode
  m_logFragmentNodeEvents = pfft_logMode
  m_countMessageNodes = 0
  m_countFragmentNodes = 0
  m_logMessageUnitNum = pfft_logUnit
  m_logFragmentUnitNum = pfft_logUnit

  !Prepare 2 lists: 1 for FLASH grid points, 1 for PFFT grid points.
  call initialise_list(pfft_flashMsgList)
  call initialise_list(pfft_pencilMsgList)


  !Create the MPI datatypes which describe the different block fragments.
  !Note: Some systems do not have an implementation for MPI_Type_create_struct()
  !which is part of MPI-2.  As such, we use MPI_Type_struct() which is safe 
  !in this particular case because we never refer to the addresses of the 
  !struct elements.  We could not (safely) get away with this on a 64 bit 
  !system if our struct consisted of different datatypes (here we would 
  !need to refer to addresses).
  offsets(1) = 0
  oldtypes(1) = FLASH_INTEGER

  !Fragment node.
  blockcounts(1) = 9 + (MDIM * 6)
  call MPI_Type_struct(1, blockCounts, offsets, oldTypes, &
       pfft_mpiMetadataType, ierr) 
  call Driver_checkMPIErrorCode(ierr)
  call MPI_Type_commit(pfft_mpiMetadataType, ierr)
  call Driver_checkMPIErrorCode(ierr)
  
end subroutine gr_pfftInitMapData
