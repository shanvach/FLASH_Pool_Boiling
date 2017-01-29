!!****if* source/Grid/GridStructures/withTriangles/pointTopoint/bitmap/gr_sbfindNeighProcs
!!
!! NAME
!!  gr_sbfindNeighProcs
!!
!! SYNOPSIS
!!
!!  gr_sbfindNeighProcs()
!!  
!! DESCRIPTION 
!!  
!!  Overview of the algoritm
!!
!!  Each processor calculates its neighboring processors and sends it to the master
!!
!!
!! ARGUMENTS 
!!
!!***

#include "constants.h"
#include "Flash.h"

  subroutine gr_sbfindNeighProcs()
    use Grid_interface, ONLY : Grid_getListOfBlocks
    use Grid_data, ONLY : gr_meshMe, gr_meshComm, gr_meshNumProcs
    use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbParticleCount, gr_sbNumBodies
    use ut_qsortInterface, ONLY : ut_qsort
    use gr_interface, ONLY : gr_findAllNeghID
    use gr_interfaceTypeDecl
    implicit none
    include "Flash_mpi.h"
    integer, allocatable, dimension(:) :: procs
    integer, allocatable, dimension(:) :: uniqueProcs
    integer, dimension(MAXBLOCKS) :: blockList
    integer :: i, j, k, b, u, p, n, ENDI, ENDJ, ENDK, allCenters, numNegh, &
         procCount, procSize, blockCount, blockID, x, ii, totalunique, jj
    type(AllBlockRegions_t) :: surrBlksSummary
    integer :: recv, send, procID, recvNeigh, sendNeigh
    integer, allocatable, dimension(:) :: neighProcs
    integer, dimension(1:2,1:gr_meshNumProcs) :: uniqueProcsCount
    integer :: sts(MPI_STATUS_SIZE,gr_meshNumProcs), recvs(gr_meshNumProcs), sends(gr_meshNumProcs), ierr, &
         recvN(gr_meshNumProcs), sendN(gr_meshNumProcs)

    recvs = MPI_REQUEST_NULL
    sends = MPI_REQUEST_NULL
    recvN = MPI_REQUEST_NULL
    sendN = MPI_REQUEST_NULL
    do x = 1, gr_sbNumBodies
       recv = 0
       send = 0
       sendNeigh = 0
       recvNeigh = 0
       if (gr_sbParticleCount(x) > 0) then
          ENDI = 3
          ENDJ = max(1,K2D*3)
          ENDK = max(1,K3D*3)
          allCenters = 2**NDIM
          call Grid_getListOfBlocks(LEAF, blockList, blockCount)
          
          !Up to "2**(NDIM-1)" neighbors for each guard cell region
          !Exactly "((ENDI*ENDJ*ENDK)-1)" guard cell regions in each block
          !Exactly "blockCount" LEAF blocks in this MPI rank.
          !Add "1" to ensure my MPI rank is in the procs array.
          procSize = (2**(NDIM-1) * ((ENDI*ENDJ*ENDK)-1) * blockCount) + 1
          allocate(procs(procSize))

          procCount = 0
          do b = 1, blockCount
             blockID = blockList(b)
             call gr_findAllNeghID(blockID, surrBlksSummary)
             do k = 1, ENDK
                do j = 1, ENDJ
                   do i = 1, ENDI
                      if ((i*j*k) /= allCenters) then
                         numNegh = surrBlksSummary % regionInfo(i,j,k) % numNegh
                         do n = 1, numNegh
                            procID = surrBlksSummary % regionInfo(i,j,k) % &
                                 details(PROCNO,n)
                            if (procID /= gr_meshMe) then
                               procCount = procCount + 1
                               procs(procCount) = procID
                            end if
                         end do
                      end if
                   end do
                end do
             end do
          end do

          if(procCount > 0) then
             call ut_qsort(procs, procCount)
             allocate(uniqueProcs(procCount+1))
             uniqueProcs = -1
             uniqueProcs(1) = gr_meshMe
             u = 1
             do p = 1, procCount
                if (.not.(any(procs(p)==uniqueProcs))) then
                   u = u + 1
                   uniqueProcs(u) = procs(p)
                end if
             end do
             deallocate(procs)
!             print *, procCount, u
!             write(*,*) "Proc",gr_meshMe,"procs", uniqueProcs
          end if
       end if
       uniqueProcsCount = 0
       totalunique = 0
       if (gr_sbBodyInfo(x) % myPE == gr_sbBodyInfo(x) % bodyMaster) then
          if (gr_sbBodyInfo(x) % sendProcs > 0) then             
             if (associated(gr_sbBodyInfo(x) % particlesPerProc)) then
                !Receive neighboring processor count and neighboring PEs for each processor
                do ii = 1,  gr_sbBodyInfo(x) % sendProcs
                   if (gr_sbBodyInfo(x) % particlesPerProc(2,ii) > 0) then
                      if (gr_sbBodyInfo(x) % particlesPerProc(1,ii) /= &
                           gr_sbBodyInfo(x) % bodyMaster) then
!                         print *, "RECV", x, gr_sbBodyInfo(x) % particlesPerProc(1,ii)
                         uniqueProcsCount(1,ii) = gr_sbBodyInfo(x) % particlesPerProc(1,ii)
                         recv = recv + 1
                         call MPI_Irecv(uniqueProcsCount(2,ii),1,FLASH_INTEGER,&
                              int(gr_sbBodyInfo(x) % particlesPerProc(1,ii)),x,gr_meshComm, recvs(recv), ierr)
                      end if
                   end if
                end do
             end if
          end if
       else
          if (gr_sbParticleCount(x) > 0) then
!             print *, "SEND", x, gr_meshMe
             send = send + 1
!             print *, "SEND", "body", x, "proc", gr_meshMe, "count", u
             call MPI_ISend(u, 1, FLASH_INTEGER, &
                  gr_sbBodyInfo(x) % bodyMaster, x, &
                  gr_meshComm, sends(send), ierr)
          end if
       end if
       if (recv > 0) then
          call MPI_WaitAll(recv, recvs, sts, ierr)
         if (ierr /= MPI_SUCCESS) then
            call Driver_abortFlash("Send MPI_Waitall error")
         endif
         do ii = 1,  gr_sbBodyInfo(x) % sendProcs
!            print *, "RECV", "body", x, "from", uniqueProcsCount(1,ii), "count", uniqueProcsCount(2,ii)
            totalunique = totalunique + uniqueProcsCount(2,ii)
         end do         
      end if

      allocate(neighProcs(totalunique))

       if (send > 0) then
          call MPI_WaitAll(send, sends, sts, ierr)
          if (ierr /= MPI_SUCCESS) then
             call Driver_abortFlash("Send MPI_Waitall error")
          endif
       end if
       
       if (gr_sbBodyInfo(x) % myPE == gr_sbBodyInfo(x) % bodyMaster) then
          if (gr_sbBodyInfo(x) % sendProcs > 0) then
             if (associated(gr_sbBodyInfo(x) % particlesPerProc)) then
                jj =1
                do ii = 1,  gr_sbBodyInfo(x) % sendProcs
                   if (gr_sbBodyInfo(x) % particlesPerProc(2,ii) > 0) then
                      if (gr_sbBodyInfo(x) % particlesPerProc(1,ii) /= &
                           gr_sbBodyInfo(x) % bodyMaster) then
                         recvNeigh = recvNeigh + 1
                         call MPI_Irecv(neighProcs(jj),&
                              uniqueProcsCount(2,ii),FLASH_INTEGER,&
                              int(gr_sbBodyInfo(x) % particlesPerProc(1,ii)),x+gr_sbNumBodies,gr_meshComm, recvN(recvNeigh), ierr)
                         jj = jj + uniqueProcsCount(2,ii)
                      end if
                   end if
                end do
                deallocate(gr_sbBodyInfo(x) % particlesPerProc)
             end if
          end if
       else
          if (gr_sbParticleCount(x) > 0) then
             sendNeigh = sendNeigh + 1
!             print *,  "Send", "body", x, "from", gr_meshMe, "procs", uniqueProcs
             call MPI_ISend(uniqueProcs, u, FLASH_INTEGER, &
                  gr_sbBodyInfo(x) % bodyMaster, x+gr_sbNumBodies, &
                  gr_meshComm, sendN(sendNeigh), ierr)
          end if
       end if

       if (sendNeigh > 0) then
          call MPI_WaitAll(sendNeigh, sendN, sts, ierr)
          if (ierr /= MPI_SUCCESS) then
             call Driver_abortFlash("Send MPI_Waitall error")
          endif          
       end if

       if (recvNeigh > 0) then
          call MPI_WaitAll(recvNeigh, recvN, sts, ierr)
          if (ierr /= MPI_SUCCESS) then
             call Driver_abortFlash("Send MPI_Waitall error")
          endif
!          print *, "recv", x, neighProcs
!          jj = 1
!          do ii = 1,  gr_sbBodyInfo(x) % sendProcs
!             print *, "RECV", "body", x, "from", uniqueProcsCount(1,ii), "procs", neighProcs(jj: uniqueProcsCount(2,ii))
!             jj = jj + uniqueProcsCount(2,ii)
!          end do

          allocate(gr_sbBodyInfo(x) % neighbors(gr_meshNumProcs))
          gr_sbBodyInfo(x) % neighbors = -1
          u = 0
          do p = 1, totalunique
             if (.not.(any(neighProcs(p)==gr_sbBodyInfo(x) % neighbors))) then
                u = u + 1
                gr_sbBodyInfo(x) % neighbors(u) = neighProcs(p)
             end if
          end do          
!          print *, x, gr_sbBodyInfo(x) % neighbors
       end if

       deallocate(neighProcs)

       if(recvNeigh > 0) then
          deallocate(gr_sbBodyInfo(x) % neighbors)
       end if

       if (gr_sbParticleCount(x) > 0) then
          deallocate(uniqueProcs)
       end if
    end do
  end subroutine gr_sbfindNeighProcs
