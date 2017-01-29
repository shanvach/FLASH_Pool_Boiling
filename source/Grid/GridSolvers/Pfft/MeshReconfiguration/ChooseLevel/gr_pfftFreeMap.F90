!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/ChooseLevel/gr_pfftFreeMap
!!
!! NAME 
!!
!! gr_pfftFreeMap
!!
!! SYNOPSIS
!!
!! gr_pfftFreeMap()
!!
!! DESCRIPTION 
!!
!! Finalises variables, arrays and data structure which were used for
!! mapping from FLASH grid -> Pencil grid and Pencil grid -> FLASH grid. 
!!
!! ARGUMENTS
!!
!! SIDE EFFECTS 
!!
!! Finalisation of module level variables, arrays, lists:
!!
!! NOTES
!! 
!! We destroy an MPI datatype in this routine which describes the metadata 
!! that we exchange between nodes.
!!
!!***

subroutine gr_pfftFreeMap()
#include "constants.h"
  use Driver_interface, ONLY : Driver_checkMPIErrorCode, &
    Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_close
  use gr_pfftMessageList, ONLY : finalise_list
  use gr_pfftData, ONLY : pfft_ndim, pfft_myPE
  use gr_pfftReconfigData, ONLY : pfft_srcGridVar, pfft_solnGridVar, &
       pfft_pfftBuf, pfft_pencilSize, pfft_procLookup, &
       pfft_mpiMetadataType, pfft_logMode, pfft_logUnit, &
       pfft_flashMsgList, pfft_pencilMsgList
  implicit none
  include "Flash_mpi.h"
  integer :: ierr, error, i

  do i = 1, pfft_ndim
     deallocate(pfft_procLookup(i) % procInfo, STAT=error)
     if (error /= 0) then
        call Driver_abortFlash("[gr_pfftReconfigFinalise]: " // &
             "Severe error. pfft_procLookup cannot be deallocated!")
     end if
  end do

  if (pfft_logMode .eqv. .true.) then
     write(pfft_logUnit,*) ""
     write(pfft_logUnit,*) "[gr_pfftFreeMap]: Freeing all metadata"
     write(pfft_logUnit,*) "-------------------------------------------------------------"
  end if

  !Destroy the metadata derived datatypes.
  call finalise_list(pfft_flashMsgList)
  call finalise_list(pfft_pencilMsgList)

  !Destroy the MPI stencils over the different fragment datatypes.
  call MPI_Type_free(pfft_mpiMetadataType, ierr)
  call Driver_checkMPIErrorCode(ierr)
  pfft_mpiMetadataType = MPI_DATATYPE_NULL

  pfft_srcGridVar = NONEXISTENT; pfft_solnGridVar = NONEXISTENT
  pfft_pencilSize = NONEXISTENT
  nullify(pfft_pfftBuf)

#ifdef LOG_PFFT
  call Logfile_close(.true.)
#endif

end subroutine gr_pfftFreeMap
