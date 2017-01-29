!!****if* source/Grid/GridStructures/withTriangles/pointTopoint/bitmap/Grid_sbSelectMaster
!!
!! NAME
!!  Grid_sbSelectMaster
!!
!! SYNOPSIS
!!
!!  Grid_sbSelectMaster()
!!  
!! DESCRIPTION 
!!  
!!  This routine is called from Grid_initDomain and from Driver_evoloveFlash
!!  for moving body. It caluclates the coordinates of each triangle for all bodies.
!!  It identifies the processor that contains the maximum no
!!  triangle centroids as the master processor.
!!
!!
!! ARGUMENTS 
!!
!!***

#include "constants.h"
#include "Flash.h"

Subroutine Grid_sbSelectMaster()
  use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkBoundBox
  use Grid_data, ONLY : gr_meshMe, gr_meshComm, gr_meshNumProcs
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, aelem, &
                        NumAelem, NumVertices, NodesPerElem
  use gr_sbInterface, only : gr_sbSendBoundBox
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Logfile_interface, ONLY : Logfile_stamp
  implicit none
  include "Flash_mpi.h"
  
  real, dimension(2,MDIM) :: boundBox
  real, dimension(MDIM) :: lowOverlap, highOverlap
  real, allocatable, dimension(:,:) :: overlapTuple, maxOverlapTuple
  integer,dimension(MAXBLOCKS) :: listOfBlocks
  integer :: count, i, jsend, jrecv,k, ierr, status(MPI_STATUS_SIZE)
  integer :: numCentroids, blkID

  logical, parameter :: logBodyStats = .true.
  integer, parameter :: MASTER = 1, SLAVE = 2, NUMSTATS=2
  integer, dimension(NUMSTATS) :: stats, minStats, maxStats, sumStats
  character(len=MAX_STRING_LENGTH) :: minStats_String, maxStats_String, &
       avgStats_String
  character(len=*), parameter, dimension(NUMSTATS) :: &
       typeStats = (/"masters/proc"," slaves/proc"/)

  integer :: OldMaster

  integer :: isrc,idest,itagv,itage, ierrv, ierre

  real, allocatable, dimension(:,:) :: DataBuf_vert
  integer, allocatable, dimension(:,:) :: Databuf_elem

  integer :: AELT_POSITION,NUM_ELEM_VARS

  real, allocatable, dimension(:,:,:) :: boundboxes
  logical :: broadcast
  integer :: sts(MPI_STATUS_SIZE,gr_meshNumProcs), recvs(gr_meshNumProcs), sends(gr_meshNumProcs)

  AELT_POSITION = AEL_POSITION  + NodesPerElem + 1
  NUM_ELEM_VARS = AELT_POSITION
  recvs = MPI_REQUEST_NULL
  sends = MPI_REQUEST_NULL

  call Timers_start("body_select_master")
  allocate(overlapTuple(2, gr_sbNumBodies))
  allocate(maxOverlapTuple(2, gr_sbNumBodies))

  call Timers_start("allreduce_boundbox")
  call Grid_getListOfBlocks(LEAF, listOfBlocks, count)

  allocate(boundboxes(LOW:HIGH, NDIM, gr_sbNumBodies))
  boundboxes=0.0
     
  do k = 1,gr_sbNumBodies
        
     if (gr_meshMe .eq. gr_sbBodyInfo(k) % bodyMaster) then
        boundboxes(LOW,1,k)  = minval(gr_sbBodyInfo(k) % xb)
        boundboxes(HIGH,1,k) = maxval(gr_sbBodyInfo(k) % xb)
#if NDIM >= 2
        boundboxes(LOW,2,k)  = minval(gr_sbBodyInfo(k) % yb)
        boundboxes(HIGH,2,k) = maxval(gr_sbBodyInfo(k) % yb)
#endif
#if NDIM == 3
        boundboxes(LOW,3,k)  = minval(gr_sbBodyInfo(k) % zb)
        boundboxes(HIGH,3,k) = maxval(gr_sbBodyInfo(k) % zb)
#endif
     endif
  enddo

  CALL MPI_ALLREDUCE(MPI_IN_PLACE, boundboxes,2*NDIM*gr_sbNumBodies,&
       FLASH_REAL, MPI_SUM,gr_meshComm,ierr)
     
  
  do k = 1,gr_sbNumBodies
     gr_sbBodyInfo(k) % boundBox(LOW,IAXIS) = boundboxes(LOW,1,k)
     gr_sbBodyInfo(k) % boundBox(HIGH,IAXIS) = boundboxes(HIGH,1,k)
#if NDIM >= 2
     gr_sbBodyInfo(k) % boundBox(LOW,JAXIS) = boundboxes(LOW,2,k)
     gr_sbBodyInfo(k) % boundBox(HIGH,JAXIS) = boundboxes(HIGH,2,k)
#endif
#if NDIM == 3
     gr_sbBodyInfo(k) % boundBox(LOW,KAXIS) = boundboxes(LOW,3,k)
     gr_sbBodyInfo(k) % boundBox(HIGH,KAXIS) = boundboxes(HIGH,3,k)
#endif        
  end do
     
  call Timers_stop("allreduce_boundbox")

  call Timers_start("find_overlap")
     
  ! Now check for overlapping space with Eulerian grid Blocks:
  do k = 1, gr_sbNumBodies
     overlapTuple(1,k) = 0.0

     do i = 1, count
        blkID = listOfBlocks(i)
        call Grid_getBlkBoundBox(blkID, boundBox)
           
        if (all(&
             gr_sbBodyInfo(k) % boundBox(LOW,1:NDIM) < boundBox(HIGH,1:NDIM) .and. &
             gr_sbBodyInfo(k) % boundBox(HIGH,1:NDIM) > boundBox(LOW,1:NDIM))) then
           lowOverlap(1:NDIM) = max(&
                gr_sbBodyInfo(k) % boundBox(LOW,1:NDIM),&
                boundBox(LOW,1:NDIM))
           
           highOverlap(1:NDIM) = min(&
                gr_sbBodyInfo(k) % boundBox(HIGH,1:NDIM),&
                boundBox(HIGH,1:NDIM))

           overlapTuple(1,k) = overlapTuple(1,k) + &
                product(highOverlap(1:NDIM) - lowOverlap(1:NDIM))
        endif
     end do

     gr_sbBodyInfo(k) % myPE = gr_meshMe
     overlapTuple(2,k) = gr_meshMe
     !The processor with blocks that overlap the most with the triangle centroids is the solid body 
     !master processor                  
  end do
  call Timers_stop("find_overlap")

  call Timers_start("load_imbalance")
  call MPI_Barrier (gr_meshComm, ierr)
  call Timers_stop("load_imbalance")

  call Timers_start("allreduce_overlap")

  call MPI_AllReduce(overlapTuple, maxOverlapTuple, gr_sbNumBodies, &
       MPI_2Double_Precision, MPI_MAXLOC, gr_meshComm, ierr)

  call Timers_stop("allreduce_overlap")
  
  call Timers_start("data_transfer")
    

  ! Now check if master changed and send body data accordingly
  do k = 1, gr_sbNumBodies

     OldMaster = gr_sbBodyInfo(k) % bodyMaster

     !write(*,*) k,'OldMaster=',OldMaster,int(maxOverlapTuple(CONSTANT_TWO,k)),gr_sbNumBodies

     gr_sbBodyInfo(k) % bodyMaster = int(maxOverlapTuple(CONSTANT_TWO,k))

     ! Did Master Processor Change? If yes send data
     if (OldMaster .ne. gr_sbBodyInfo(k) % bodyMaster) then

        NumVertices  = gr_sbBodyInfo(k)%NumVertices 
        NumAelem     = gr_sbBodyInfo(k)%NumAelem
        jsend = 0
        jrecv = 0

        ! Communication info
        isrc  = OldMaster
        idest = gr_sbBodyInfo(k)%bodyMaster
        itagv =  k
        itage =  k + gr_sbNumBodies

        if(gr_meshMe .eq. gr_sbBodyInfo(k) % bodyMaster) then  ! I'm the new bodyMaster
           jrecv = jrecv + 1

           ! Then allocate the other fields in gr_sbBodyInfo
           allocate(gr_sbBodyInfo(k) % xbus(NumVertices), gr_sbBodyInfo(k) % ybus(NumVertices))
           allocate(gr_sbBodyInfo(k) % xb(NumVertices), gr_sbBodyInfo(k) % yb(NumVertices)) !Vertex points
           allocate(gr_sbBodyInfo(k) % ubd(NumVertices), gr_sbBodyInfo(k) % vbd(NumVertices)) !Vertex velocities
           allocate(gr_sbBodyInfo(k) % ubdd(NumVertices), gr_sbBodyInfo(k) % vbdd(NumVertices)) !Vertex accelrations
           allocate(gr_sbBodyInfo(k) % AAREANODE(NumVertices), gr_sbBodyInfo(k) % AANGNODE(NumVertices))
           allocate(gr_sbBodyInfo(k) % dutdn(NumVertices), gr_sbBodyInfo(k) % dvtdn(NumVertices))
           allocate(gr_sbBodyInfo(k) % nxL(NumVertices), gr_sbBodyInfo(k) % nyL(NumVertices))   !  nzL for 3D
           allocate(gr_sbBodyInfo(k) % AAREAELEM(NumAelem))
           allocate(gr_sbBodyinfo(k) % AELEMNODE(NumVertices,10))    ! Array of elements connected to each node 
           allocate(gr_sbBodyInfo(k) % AELEM(NodesPerElem + CONSTANT_ONE,NumAelem))        ! Array of connectivity for each body
           allocate(gr_sbBodyInfo(k) % sbPtNumXi(NumAelem),gr_sbBodyInfo(k) % sbPtNumEta(NumAelem))
           allocate(gr_sbBodyInfo(k) % sbPtNumElem(NumAelem))
           allocate(gr_sbBodyInfo(k) % AELTYPE(NumAelem))
#if NDIM == 2
           allocate(gr_sbBodyInfo(k) % sb(NumVertices))
#endif

#if NDIM == 3
           allocate(gr_sbBodyInfo(k) % zbus(NumVertices))
           allocate(gr_sbBodyInfo(k) % zb(NumVertices))
           allocate(gr_sbBodyInfo(k) % wbd(NumVertices))  ! vertex velocity
           allocate(gr_sbBodyInfo(k) % wbdd(NumVertices)) ! vertex acceleration
           allocate(gr_sbBodyInfo(k) % dwtdn(NumVertices))
           allocate(gr_sbBodyInfo(k) % nzL(NumAelem))
#endif

           ! Now Receive beast buffers:
           allocate(DataBuf_vert(NumVertices,NUM_VERT_VARS))
           allocate(Databuf_elem(NUM_ELEM_VARS,NumAelem))
           !DataBuf_vert(:,:) = 0.
           !DataBuf_elem(:,:) = 0.

           call mpi_Irecv(DataBuf_vert , NumVertices*NUM_VERT_VARS, FLASH_REAL, isrc, itagv,gr_meshComm,recvs(jrecv),ierrv)
           call mpi_Irecv(Databuf_elem , NumAelem*NUM_ELEM_VARS, FLASH_INTEGER, isrc, itage,gr_meshComm,recvs(jrecv),ierre)

        elseif (gr_meshMe .eq. OldMaster) then ! I send the stuff to the new Master
           jsend = jsend + 1
           ! Create beast buffers:
           allocate(DataBuf_vert(NumVertices,NUM_VERT_VARS))
           allocate(Databuf_elem(NUM_ELEM_VARS,NumAelem))           
           DataBuf_vert(:,:) = 0.
           DataBuf_elem(:,:) = 0.
!           gr_sbBodyInfo(k) % xbus(1) = 1.0
!           gr_sbBodyInfo(k) % xbus(2) = 2.0
           ! Vertex related data:
           DataBuf_vert(1:NumVertices,XBUS_POSITION)  =  gr_sbBodyInfo(k) % xbus(1:NumVertices)
           DataBuf_vert(1:NumVertices,YBUS_POSITION)  =  gr_sbBodyInfo(k) % ybus(1:NumVertices)
           DataBuf_vert(1:NumVertices,XB_POSITION)    =  gr_sbBodyInfo(k) % xb(1:NumVertices)  
           DataBuf_vert(1:NumVertices,YB_POSITION)    =  gr_sbBodyInfo(k) % yb(1:NumVertices)  
           DataBuf_vert(1:NumVertices,UBD_POSITION)   =  gr_sbBodyInfo(k) % ubd(1:NumVertices) 
           DataBuf_vert(1:NumVertices,VBD_POSITION)   =  gr_sbBodyInfo(k) % vbd(1:NumVertices) 
           DataBuf_vert(1:NumVertices,UBDD_POSITION)  =  gr_sbBodyInfo(k) % ubdd(1:NumVertices)
           DataBuf_vert(1:NumVertices,VBDD_POSITION)  =  gr_sbBodyInfo(k) % vbdd(1:NumVertices)
           DataBuf_vert(1:NumVertices,NXL_POSITION)   =  gr_sbBodyInfo(k) % nxL(1:NumVertices) 
           DataBuf_vert(1:NumVertices,NYL_POSITION)   =  gr_sbBodyInfo(k) % nyL(1:NumVertices) 
#if NDIM == 2                                                                                          
           DataBuf_vert(1:NumVertices,SB_POSITION)    =  gr_sbBodyInfo(k) % sb(1:NumVertices)  
#elif NDIM == 3                                                                                         
           DataBuf_vert(1:NumVertices,ZBUS_POSITION)  =  gr_sbBodyInfo(k) % zbus(1:NumVertices)
           DataBuf_vert(1:NumVertices,ZB_POSITION)    =  gr_sbBodyInfo(k) % zb(1:NumVertices)  
           DataBuf_vert(1:NumVertices,WBD_POSITION)   =  gr_sbBodyInfo(k) % wbd(1:NumVertices) 
           DataBuf_vert(1:NumVertices,WBDD_POSITION)  =  gr_sbBodyInfo(k) % wbdd(1:NumVertices)
           DataBuf_vert(1:NumVertices,NZL_POSITION)   =  gr_sbBodyInfo(k) % nzL(1:NumVertices) 
#endif

           ! Send data to new Master
 !          print *, "body", k, "old master", gr_sbBodyInfo(k) % xbus(1), gr_sbBodyInfo(k) % xbus(2)
           call mpi_Isend(DataBuf_vert , NumVertices*NUM_VERT_VARS, FLASH_REAL, idest, itagv,gr_meshComm,sends(jsend),ierrv)

           ! Elem related data:
           Databuf_elem(AEL_POSITION:AEL_POSITION+NodesPerElem,1:NumAelem)  = gr_sbBodyInfo(k) % AELEM(1:NodesPerElem+CONSTANT_ONE,1:NumAelem)
           Databuf_elem(AELT_POSITION,1:NumAelem) = gr_sbBodyInfo(k) % AELTYPE(1:NumAelem)

           ! Send data to new Master
           call mpi_Isend(Databuf_elem , NumAelem*NUM_ELEM_VARS, FLASH_INTEGER, idest, itage,gr_meshComm,sends(jsend),ierre)
        endif

        if (jrecv > 0) then
           call MPI_WaitAll(jrecv, recvs, sts, ierr)
           if (ierr /= MPI_SUCCESS) then
              call Driver_abortFlash("Send MPI_Waitall error")
           endif
           gr_sbBodyInfo(k) % xbus(1:NumVertices) = DataBuf_vert(1:NumVertices,XBUS_POSITION)
!           print *, "body", k, "new master", DataBuf_vert(1,XBUS_POSITION), DataBuf_vert(2,XBUS_POSITION)
           gr_sbBodyInfo(k) % ybus(1:NumVertices) = DataBuf_vert(1:NumVertices,YBUS_POSITION)
           gr_sbBodyInfo(k) % xb(1:NumVertices)   = DataBuf_vert(1:NumVertices,XB_POSITION)
           gr_sbBodyInfo(k) % yb(1:NumVertices)   = DataBuf_vert(1:NumVertices,YB_POSITION)
           gr_sbBodyInfo(k) % ubd(1:NumVertices)  = DataBuf_vert(1:NumVertices,UBD_POSITION)
           gr_sbBodyInfo(k) % vbd(1:NumVertices)  = DataBuf_vert(1:NumVertices,VBD_POSITION)
           gr_sbBodyInfo(k) % ubdd(1:NumVertices) = DataBuf_vert(1:NumVertices,UBDD_POSITION)
           gr_sbBodyInfo(k) % vbdd(1:NumVertices) = DataBuf_vert(1:NumVertices,VBDD_POSITION)
           gr_sbBodyInfo(k) % nxL(1:NumVertices)  = DataBuf_vert(1:NumVertices,NXL_POSITION)
           gr_sbBodyInfo(k) % nyL(1:NumVertices)  = DataBuf_vert(1:NumVertices,NYL_POSITION)
#if NDIM == 2
           gr_sbBodyInfo(k) % sb(1:NumVertices)   = DataBuf_vert(1:NumVertices,SB_POSITION)
#elif NDIM == 3
           gr_sbBodyInfo(k) % zbus(1:NumVertices) = DataBuf_vert(1:NumVertices,ZBUS_POSITION)
           gr_sbBodyInfo(k) % zb(1:NumVertices)   = DataBuf_vert(1:NumVertices,ZB_POSITION)
           gr_sbBodyInfo(k) % wbd(1:NumVertices)  = DataBuf_vert(1:NumVertices,WBD_POSITION)
           gr_sbBodyInfo(k) % wbdd(1:NumVertices) = DataBuf_vert(1:NumVertices,WBDD_POSITION)
           gr_sbBodyInfo(k) % nzL(1:NumVertices)  = DataBuf_vert(1:NumVertices,NZL_POSITION)
#endif
           gr_sbBodyInfo(k) % AELEM(1:NodesPerElem+CONSTANT_ONE,1:NumAelem) = Databuf_elem(AEL_POSITION:AEL_POSITION+NodesPerElem,1:NumAelem)
           gr_sbBodyInfo(k) % AELTYPE(1:NumAelem) = Databuf_elem(AELT_POSITION,1:NumAelem)
           deallocate(Databuf_vert)
           deallocate(DataBuf_elem)
        end if

         if (jsend > 0) then
           call MPI_WaitAll(jsend, sends, sts, ierr)
           if (ierr /= MPI_SUCCESS) then
              call Driver_abortFlash("Send MPI_Waitall error")
           endif
           deallocate(Databuf_vert)
           deallocate(DataBuf_elem)
           deallocate(gr_sbBodyInfo(k) % xbus, gr_sbBodyInfo(k) % ybus)
           deallocate(gr_sbBodyInfo(k) % xb, gr_sbBodyInfo(k) % yb) !Vertex points 
           deallocate(gr_sbBodyInfo(k) % ubd, gr_sbBodyInfo(k) % vbd) !Vertex velocities 
           deallocate(gr_sbBodyInfo(k) % ubdd, gr_sbBodyInfo(k) % vbdd) !Vertex accelrations
           deallocate(gr_sbBodyInfo(k) % AAREANODE, gr_sbBodyInfo(k) % AANGNODE)
           deallocate(gr_sbBodyInfo(k) % dutdn, gr_sbBodyInfo(k) % dvtdn)
           deallocate(gr_sbBodyInfo(k) % nxL, gr_sbBodyInfo(k) % nyL)   !  nzL for 3D 
           deallocate(gr_sbBodyInfo(k) % AAREAELEM)
           deallocate(gr_sbBodyinfo(k) % AELEMNODE)    ! Array of elements connected to each node 
           deallocate(gr_sbBodyInfo(k) % AELEM)        ! Array of connectivity for each body 
           deallocate(gr_sbBodyInfo(k) % sbPtNumXi,gr_sbBodyInfo(k) % sbPtNumEta)
           deallocate(gr_sbBodyInfo(k) % sbPtNumElem)
           deallocate(gr_sbBodyInfo(k) % AELTYPE)
#if NDIM == 2
           deallocate(gr_sbBodyInfo(k) % sb)
#endif

#if NDIM == 3
           deallocate(gr_sbBodyInfo(k) % zbus)
           deallocate(gr_sbBodyInfo(k) % zb)
           deallocate(gr_sbBodyInfo(k) % wbd)  ! vertex velocity 
           deallocate(gr_sbBodyInfo(k) % wbdd) ! vertex accelaration 
           deallocate(gr_sbBodyInfo(k) % dwtdn)
           deallocate(gr_sbBodyInfo(k) % nzL)
#endif
        end if
     endif

!     call MPI_Barrier (gr_meshComm, ierr)

  enddo
  call Timers_stop("data_transfer")

  ! Do log stats for Bodies:
  stats(:) = 0
  do k = 1, gr_sbNumBodies
     if (gr_sbBodyInfo(k) % bodyMaster == gr_sbBodyInfo(k) % myPE) then
        stats(MASTER) = stats(MASTER) + 1
     else
        if (overlapTuple(1,k) > 0.0) then
           stats(SLAVE) = stats(SLAVE) + 1
        end if
     end if
  enddo


  if (logBodyStats) then
     call MPI_AllReduce(stats, minStats, NUMSTATS, FLASH_INTEGER, MPI_MIN, &
          gr_meshComm, ierr)
     call MPI_AllReduce(stats, maxStats, NUMSTATS, FLASH_INTEGER, MPI_MAX, &
          gr_meshComm, ierr)
     call MPI_AllReduce(stats, sumStats, NUMSTATS, FLASH_INTEGER, MPI_SUM, &
          gr_meshComm, ierr)

     do i = 1, NUMSTATS
        if (gr_meshMe == 0) then
           write (minStats_string, '(i10)') minStats(i)
           write (maxStats_string, '(i10)') maxStats(i)
           write (avgStats_string, '(i10)') int(sumStats(i) / real(gr_meshNumProcs))
           call Logfile_stamp(typeStats(i) // & 
                trim(minStats_string) // ' (min)  ' // &
                trim(maxStats_string) // ' (max)  ' // &
                trim(avgStats_string) // ' (avg) ', tag="body info")
        end if
     end do
  end if

  deallocate(overlapTuple)
  deallocate(maxOverlapTuple)
  deallocate(boundboxes)
  call Timers_stop("body_select_master")

End Subroutine Grid_sbSelectMaster
