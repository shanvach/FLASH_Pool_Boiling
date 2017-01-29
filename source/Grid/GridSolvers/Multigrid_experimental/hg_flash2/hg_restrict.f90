!******************************************************************************

! Routine:      hg_restrict

! Description:  Restrict data from blocks on a given level to their parents,
!               for use in constructing the composed solution on the parents.
!               The parent blocks are first zeroed, and then the restrictions
!               are performed.

! Parameters:   ifrom        Variable index for variable to restrict from
!               ito          Variable index for variable to restrict to
!               level        Level from which to restrict


subroutine hg_restrict(level, ifrom, ito)

!==============================================================================

use perfmon
use dBase, ONLY:  dBaseRefinementLevel, dBaseGetDataPtrSingleBlock, &
                  dBaseReleaseDataPtrSingleBlock, GC, nxb, nyb, nzb, &
                  dBaseBlockSize, ndim, nguard, dBasePropertyInteger,&
                  dBaseNodeType, maxblocks, nchild, k2d, k3d
use dBaseDeclarations, ONLY: parent, child

implicit none

include 'mpif.h'

integer, intent(in)          :: ifrom, ito, level

integer                      :: b, c, p, h, i, j, k, lnblocks
integer                      :: ierr1, ierr2, ierr3
integer                      :: i1, i2, j1, j2, k1, k2, ichild
integer                      :: ierr, nsent
real, pointer                :: solnData(:,:,:,:)
integer                      :: status(MPI_STATUS_SIZE)

integer, allocatable, save   :: send_req(:)
real, allocatable, save      :: send_parent_data(:,:,:,:)
real, allocatable, save      :: recv_parent_data(:,:,:)

integer, save                :: n1, n2, n3, mype, nbbuf
logical, save                :: first_call = .true.

!==============================================================================

call timer_start("hg_restrict")

if (first_call) then
  n1 = max(nxb/2, 1)
  n2 = max(nyb/2, 1)
  n3 = max(nzb/2, 1)
  nbbuf = max(maxblocks, 1)
  allocate(send_parent_data(n1,n2,n3,nbbuf), stat=ierr1)
  allocate(recv_parent_data(n1,n2,n3), stat=ierr2)
  allocate(send_req(nbbuf), stat=ierr3)
  if ((ierr1 /= 0) .or. (ierr2 /= 0) .or. (ierr3 /= 0)) &
    call abort_flash("Allocation error in hg_restrict")
  mype = dBasePropertyInteger("MyProcessor")
  first_call = .false.
endif

lnblocks = dBasePropertyInteger("LocalNumberOfBlocks")

nsent = 0
h = 1

do b = 1, lnblocks
  if (dBaseRefinementLevel(b) == level) then

    solnData => dBaseGetDataPtrSingleBlock(b, GC)

! First, compute restricted values for this block.

    if (ndim == 1) then

      do i = 1, n1
        i1 = nguard + 2*i - 1
        i2 = i1 + 1
        send_parent_data(i,1,1,h)  = 0.5 * (solnData(ifrom,i1,1,1) + &
                                            solnData(ifrom,i2,1,1))
      enddo

    else if (ndim == 2) then

      do j = 1, n2
        j1 = nguard + 2*j - 1
        j2 = j1 + 1
        do i = 1, n1
          i1 = nguard + 2*i - 1
          i2 = i1 + 1
          send_parent_data(i,j,1,h) = 0.25*(solnData(ifrom,i1,j1,1)+&
                                            solnData(ifrom,i2,j1,1)+&
                                            solnData(ifrom,i1,j2,1)+&
                                            solnData(ifrom,i2,j2,1))
        enddo
      enddo

    else ! ndim == 3

      do k = 1, n3
        k1 = nguard + 2*k - 1
        k2 = k1 + 1
        do j = 1, n2
          j1 = nguard + 2*j - 1
          j2 = j1 + 1
          do i = 1, n1
            i1 = nguard + 2*i - 1
            i2 = i1 + 1
            send_parent_data(i,j,k,h) = 0.125*(solnData(ifrom,i1,j1,k1)+&
                                               solnData(ifrom,i2,j1,k1)+&
                                               solnData(ifrom,i1,j2,k1)+&
                                               solnData(ifrom,i2,j2,k1)+&
                                               solnData(ifrom,i1,j1,k2)+&
                                               solnData(ifrom,i2,j1,k2)+&
                                               solnData(ifrom,i1,j2,k2)+&
                                               solnData(ifrom,i2,j2,k2))
          enddo
        enddo
      enddo

    endif

    call dBaseReleaseDataPtrSingleBlock(b, solnData)

! Next, if parent is on this processor, copy the restricted values directly.
! If parent is off-processor, send the data to the owning processor
! (non-blocking).

    if (parent(2,b) == mype) then ! local parent
      p = parent(1,b)
      do ichild = 1, nchild
        if (child(1,ichild,p) == b) then
          c = ichild
          exit
        endif
      enddo
      if ((c == 1) .or. (c == 3) .or. (c == 5) .or. (c == 7)) then
        i1 = nguard + 1
        i2 = nguard + n1
      else
        i1 = nguard + nxb - n1 + 1
        i2 = nguard + nxb
      endif
      if ((c == 1) .or. (c == 2) .or. (c == 5) .or. (c == 6)) then
        j1 = nguard*k2d + 1
        j2 = nguard*k2d + n2
      else
        j1 = nguard*k2d + nyb - n2 + 1
        j2 = nguard*k2d + nyb
      endif
      if ((c == 1) .or. (c == 2) .or. (c == 3) .or. (c == 4)) then
        k1 = nguard*k3d + 1
        k2 = nguard*k3d + n3
      else
        k1 = nguard*k3d + nzb - n3 + 1
        k2 = nguard*k3d + nzb
      endif
      solnData => dBaseGetDataPtrSingleBlock(p, GC)
      solnData(ito,i1:i2,j1:j2,k1:k2) = send_parent_data(:,:,:,h)
      call dBaseReleaseDataPtrSingleBlock(p, solnData)
    else                          ! remote parent
      nsent = nsent + 1
      call mpi_isend(send_parent_data(1,1,1,h), n1*n2*n3, &
                     MPI_DOUBLE_PRECISION, parent(2,b), b, &
                     MPI_COMM_WORLD, send_req(nsent), ierr)
      h = h + 1
    endif

    if (h > nbbuf) &
      call abort_flash("Buffer space exceeded in hg_restrict")

  endif
enddo

do b = 1, lnblocks
  if ((dBaseRefinementLevel(b) == level-1) .and. (dBaseNodeType(b) > 1)) then
    do c = 1, nchild
      if (child(2,c,b) /= mype) then

! If child is on another processor, receive the restricted data from
! the child (blocking).

        call mpi_recv(recv_parent_data(1,1,1), n1*n2*n3, &
                      MPI_DOUBLE_PRECISION, child(2,c,b), child(1,c,b), &
                      MPI_COMM_WORLD, status, ierr)

        if ((c == 1) .or. (c == 3) .or. (c == 5) .or. (c == 7)) then
          i1 = nguard + 1
          i2 = nguard + n1
        else
          i1 = nguard + nxb - n1 + 1
          i2 = nguard + nxb
        endif
        if ((c == 1) .or. (c == 2) .or. (c == 5) .or. (c == 6)) then
          j1 = nguard*k2d + 1
          j2 = nguard*k2d + n2
        else
          j1 = nguard*k2d + nyb - n2 + 1
          j2 = nguard*k2d + nyb
        endif
        if ((c == 1) .or. (c == 2) .or. (c == 3) .or. (c == 4)) then
          k1 = nguard*k3d + 1
          k2 = nguard*k3d + n3
        else
          k1 = nguard*k3d + nzb - n3 + 1
          k2 = nguard*k3d + nzb
        endif

        solnData => dBaseGetDataPtrSingleBlock(b, GC)
        solnData(ito,i1:i2,j1:j2,k1:k2) = recv_parent_data(:,:,:)
        call dBaseReleaseDataPtrSingleBlock(b, solnData)

      endif
    enddo
  endif
enddo

call mpi_waitall(nsent, send_req, MPI_STATUSES_IGNORE, ierr)

call timer_stop("hg_restrict")

!==============================================================================

return
end subroutine hg_restrict
