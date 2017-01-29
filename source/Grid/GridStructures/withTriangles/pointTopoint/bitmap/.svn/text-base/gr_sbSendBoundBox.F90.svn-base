!!****if* source/Grid/GridStructures/withTriangles/pointTopoint/bitmap/gr_sbSendBoundBox
!!
!! NAME
!!  gr_sbSendBoundBox
!!
!! SYNOPSIS
!!
!!  call gr_sbSendBoundBox
!!
!! DESCRIPTION
!!
!!  Called from gr_sbInit. Calculate the boundbox of each solid body.  
!!
!! ARGUMENTS

#include "constants.h"
#include "Flash.h"

Subroutine gr_sbSendBoundBox()
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies
  use Grid_data, ONLY : gr_meshComm, gr_meshNumProcs, gr_meshMe
 
  implicit none
  include "Flash_mpi.h"
  integer :: b, iSend, iRecv, ierr, i
  integer, dimension(gr_meshNumProcs) :: sendreq
  integer, dimension(MPI_STATUS_SIZE, gr_meshNumProcs) :: sendstatus
  integer, dimension(gr_meshNumProcs) :: recvreq
  integer, dimension(MPI_STATUS_SIZE, gr_meshNumProcs) :: recvstatus

  sendreq = MPI_REQUEST_NULL
  recvreq = MPI_REQUEST_NULL
  do b = 1, gr_sbNumBodies   
     iRecv = 0
     iSend = 0
     if (gr_meshMe == gr_sbBodyInfo(b) % bodyMaster) then
        do i = 0, gr_meshNumProcs-1
           if (i /= gr_sbBodyInfo(b) % bodyMaster) then
              iSend = iSend + 1
              call MPI_Isend(gr_sbBodyInfo(b) % boundBox(:,:), MDIM*2,&
                   FLASH_REAL,i,b,gr_meshComm,sendreq(iSend),ierr)
           endif
        enddo
     else
        iRecv = iRecv + 1
        call MPI_Irecv(gr_sbBodyInfo(b) % boundBox(:,:), MDIM*2,&
             FLASH_REAL, gr_sbBodyInfo(b) % bodyMaster,b,gr_meshComm,&
             recvreq(iRecv),ierr)
     endif
     if (iSend > 0) then
        call MPI_WAITALL(iSend, sendreq, sendstatus, ierr)
        if (ierr /= MPI_SUCCESS) then
           call Driver_abortFlash("Send MPI_Waitall error")
        endif
     endif

     if (iRecv > 0) then
        call MPI_WAITALL(iRecv, recvreq, recvstatus, ierr)
        if (ierr /= MPI_SUCCESS) then
           call Driver_abortFlash("Recv MPI_Waitall error")
        endif
     endif
  enddo
!  return
End Subroutine gr_sbSendBoundBox
