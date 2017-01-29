!******************************************************************************

! Routine:      hg_prolong_bndries

! Description:  Prolongate (interpolate) boundary information from parent
!               blocks on a given level of refinement to their children at
!               the next finer level.

!               For any child block we generally have two types of boundary:
!               those coinciding with one of the parent block's boundaries,
!               and those passing through the interior of the parent block.
!               For the first type of boundary, we interpolate surface
!               boundary values taken from the parent's first layer of solution
!               variable guard cell entries.  For the second type of boundary,
!               we average the cell averages straddling the boundary in the
!               parent block to estimate surface boundary values, then
!               interpolate these estimates to the child.  All boundary values
!               are written to the first layer of guard cell entries on the
!               child block.

! Parameters:   ifrom        Variable index for variable to interpolate from
!               ito          Variable index for variable to interpolate to
!               level        Level from which to interpolate
!               ichild       If /= 0, only prolongate boundary data for the
!                            indicated child index.


subroutine hg_prolong_bndries(level, ifrom, ito, ichild)

!==============================================================================

use dBase, ONLY:  dBaseRefinementLevel, dBaseGetDataPtrSingleBlock, &
                  dBaseReleaseDataPtrSingleBlock, GC, nxb, nyb, nzb, &
                  dBaseBlockSize, ndim, nguard, k2d, k3d, nchild, nfaces, &
                  dBasePropertyInteger, dBaseNodeType, maxblocks, &
                  dBaseNeighborBlockList
use dBaseDeclarations, ONLY: child, parent, work

use mg_common
use perfmon

implicit none

include 'mpif.h'

integer, intent(in)          :: ifrom, ito, level, ichild

integer                      :: b, c, h, i, j, k, lnblocks, ierr1, ierr2, ierr3
integer                      :: i1, i2, j1, j2, ii, jj, kk
integer                      :: ierr, nsent
logical                      :: any_sent
real, pointer                :: solnData(:,:,:,:)
integer                      :: status(MPI_STATUS_SIZE)

integer, allocatable, save   :: send_req(:)
real, allocatable, save      :: send_child_data(:,:,:,:,:)
real, allocatable, save      :: recv_child_data(:,:,:)

integer, save                :: nmax1, nmax2, mype, nbbuf
integer, save                :: n1off(8), n2off(8), n3off(8)
logical, save                :: first_call = .true.

real, save                   :: Px(-2:2,nxb,nchild), Py(-2:2,nyb,nchild), Pz(-2:2,nzb,nchild)
real, save                   :: Pns(-2:2,2), Pew(-2:2,2), Pud(-2:2,2)

integer                      :: nbr_blks(6)

!==============================================================================

call timer_start("hg_prolong_bndries")

if (first_call) then
  if (ndim == 1) then
    nmax1 = 1
    nmax2 = 1
  elseif (ndim == 2) then
    nmax1 = max(nxb, nyb)
    nmax2 = 1
  else ! ndim == 3
    nmax1 = max(nxb, nyb)
    nmax2 = max(nyb, nzb)
  endif
  nbbuf = max(maxblocks/8, 1)
  allocate(send_child_data(nmax1,nmax2,nfaces,nchild,nbbuf), stat=ierr1)
  allocate(recv_child_data(nmax1,nmax2,nfaces), stat=ierr2)
  allocate(send_req(nchild*nbbuf), stat=ierr3)
  if ((ierr1 /= 0) .or. (ierr2 /= 0) .or. (ierr3 /= 0)) &
    call abort_flash("Allocation error in hg_prolong_bndries")
  mype = dBasePropertyInteger("MyProcessor")
  n1off(1) = 0
  n2off(1) = 0
  n3off(1) = 0
  n1off(2) = nxb/2
  n2off(2) = 0
  n3off(2) = 0
  n1off(3) = 0
  n2off(3) = nyb/2
  n3off(3) = 0
  n1off(4) = nxb/2
  n2off(4) = nyb/2
  n3off(4) = 0
  n1off(5) = 0
  n2off(5) = 0
  n3off(5) = nzb/2
  n1off(6) = nxb/2
  n2off(6) = 0
  n3off(6) = nzb/2
  n1off(7) = 0
  n2off(7) = nyb/2
  n3off(7) = nzb/2
  n1off(8) = nxb/2
  n2off(8) = nyb/2
  n3off(8) = nzb/2

  if (ndim == 1) then

    Pew(:,1) = (/ -1./12.,  7./12.,  7./12., -1./12.,  0.    /)
    Pew(:,2) = (/  0.,     -1./12.,  7./12.,  7./12., -1./12. /)

  endif

  if (ndim == 2) then

    do c = 1, 4
      do i = 1, nxb-1, 2
        Px(:,i,  c) = (/ -3./128.,  11./64., 1., -11./64.,  3./128. /)
        Px(:,i+1,c) = (/  3./128., -11./64., 1.,  11./64., -3./128. /)
      enddo
      do i = 1, nyb-1, 2
        Py(:,i,  c) = (/ -3./128.,  11./64., 1., -11./64.,  3./128. /)
        Py(:,i+1,c) = (/  3./128., -11./64., 1.,  11./64., -3./128. /)
      enddo
    enddo

    Pns(:,1) = (/ -1./12.,  7./12.,  7./12., -1./12.,  0.    /)
    Pns(:,2) = (/  0.,     -1./12.,  7./12.,  7./12., -1./12. /)
    Pew(:,1) = (/ -1./12.,  7./12.,  7./12., -1./12.,  0.    /)
    Pew(:,2) = (/  0.,     -1./12.,  7./12.,  7./12., -1./12. /)

  endif

  if (ndim == 3) then

    do c = 1, 8
      do i = 1, nxb-1, 2
        Px(:,i,  c) = (/ -3./128.,  11./64., 1., -11./64.,  3./128. /)
        Px(:,i+1,c) = (/  3./128., -11./64., 1.,  11./64., -3./128. /)
      enddo
      do i = 1, nyb-1, 2
        Py(:,i,  c) = (/ -3./128.,  11./64., 1., -11./64.,  3./128. /)
        Py(:,i+1,c) = (/  3./128., -11./64., 1.,  11./64., -3./128. /)
      enddo
      do i = 1, nzb-1, 2
        Pz(:,i,  c) = (/ -3./128.,  11./64., 1., -11./64.,  3./128. /)
        Pz(:,i+1,c) = (/  3./128., -11./64., 1.,  11./64., -3./128. /)
      enddo
    enddo

    Pns(:,1) = (/ -1./12.,  7./12.,  7./12., -1./12.,  0.    /)
    Pns(:,2) = (/  0.,     -1./12.,  7./12.,  7./12., -1./12. /)
    Pew(:,1) = (/ -1./12.,  7./12.,  7./12., -1./12.,  0.    /)
    Pew(:,2) = (/  0.,     -1./12.,  7./12.,  7./12., -1./12. /)
    Pud(:,1) = (/ -1./12.,  7./12.,  7./12., -1./12.,  0.    /)
    Pud(:,2) = (/  0.,     -1./12.,  7./12.,  7./12., -1./12. /)

  endif

  first_call = .false.
endif

! Use (EXCHANGE_WORK, CONTINUE_SERIES) under the assumption that we're
! calling hg_prolong_bndries immediately after a call to hg_solve_level
! or hg_residual -- performance savings by not filling work() again

if (ifrom == mg_soln_index) then
  call mg_bndry(level, ifrom, nguard, 0, EXCHANGE_WORK, STANDALONE, 1)
else
  call mg_bndry(level, ifrom, nguard, 0, EXCHANGE_WORK, STANDALONE, 0)
endif

lnblocks = dBasePropertyInteger("LocalNumberOfBlocks")

nsent = 0
h = 1

do b = 1, lnblocks

  if ((dBaseRefinementLevel(b) == level) .and. (dBaseNodeType(b) > 1)) then

    solnData => dBaseGetDataPtrSingleBlock(b, GC)

    nbr_blks = dBaseNeighborBlockList(b)

! First, compute interpolated boundary values for each of this block's
! children.

    if (ndim == 1) then

      do c = 1, nchild
        if ((ichild == 0) .or. (ichild == c)) then

          send_child_data(:,:,:,c,h) = 0.

          do ii = -2, 2
            send_child_data(1,1,1,c,h) = send_child_data(1,1,1,c,h) + Pew(ii,1)*work(nguard+n1off(c)+1+ii,1,1,b)
            send_child_data(1,1,2,c,h) = send_child_data(1,1,2,c,h) + Pew(ii,2)*work(nguard+n1off(c)+nxb/2+ii,1,1,b)
          enddo

        endif
      enddo

    else if (ndim == 2) then

      do c = 1, nchild
        if ((ichild == 0) .or. (ichild == c)) then

          send_child_data(:,:,:,c,h) = 0.

          ! N & S edges

          do i = 1, nxb
            do ii = -2, 2
              do jj = -2, 2
                send_child_data(i,1,3,c,h) = send_child_data(i,1,3,c,h) + &
                    Px(ii,i,c)*Pns(jj,1)*work(nguard+n1off(c)+1+(i-1)/2+ii,nguard+n2off(c)+1+jj,1,b)
                send_child_data(i,1,4,c,h) = send_child_data(i,1,4,c,h) + &
                    Px(ii,i,c)*Pns(jj,2)*work(nguard+n1off(c)+1+(i-1)/2+ii,nguard+n2off(c)+nyb/2+jj,1,b)
              enddo
            enddo
          enddo

          ! E & W edges

          do j = 1, nyb
            do jj = -2, 2
              do ii = -2, 2
                send_child_data(j,1,1,c,h) = send_child_data(j,1,1,c,h) + &
                    Py(jj,j,c)*Pew(ii,1)*work(nguard+n1off(c)+1+ii,nguard+n2off(c)+1+(j-1)/2+jj,1,b)
                send_child_data(j,1,2,c,h) = send_child_data(j,1,2,c,h) + &
                    Py(jj,j,c)*Pew(ii,2)*work(nguard+n1off(c)+nxb/2+ii,nguard+n2off(c)+1+(j-1)/2+jj,1,b)
              enddo
            enddo
          enddo

        endif
      enddo

    else ! ndim == 3

      do c = 1, nchild
        if ((ichild == 0) .or. (ichild == c)) then

          send_child_data(:,:,:,c,h) = 0.

          ! N & S faces

          do k = 1, nzb
            do i = 1, nxb
              do kk = -2, 2
                do ii = -2, 2
                  do jj = -2, 2
                    send_child_data(i,k,3,c,h) = send_child_data(i,k,3,c,h) + &
                        Px(ii,i,c)*Pz(kk,k,c)*Pns(jj,1)*work(nguard+n1off(c)+1+(i-1)/2+ii,nguard+n2off(c)+1+jj,nguard+n3off(c)+1+(k-1)/2+kk,b)
                    send_child_data(i,k,4,c,h) = send_child_data(i,k,4,c,h) + &
                        Px(ii,i,c)*Pz(kk,k,c)*Pns(jj,2)*work(nguard+n1off(c)+1+(i-1)/2+ii,nguard+n2off(c)+nyb/2+jj,nguard+n3off(c)+1+(k-1)/2+kk,b)
                  enddo
                enddo
              enddo
            enddo
          enddo

          ! E & W faces

          do k = 1, nzb
            do j = 1, nyb
              do kk = -2, 2
                do jj = -2, 2
                  do ii = -2, 2
                    send_child_data(j,k,1,c,h) = send_child_data(j,k,1,c,h) + &
                        Py(jj,j,c)*Pz(kk,k,c)*Pew(ii,1)*work(nguard+n1off(c)+1+ii,nguard+n2off(c)+1+(j-1)/2+jj,nguard+n3off(c)+1+(k-1)/2+kk,b)
                    send_child_data(j,k,2,c,h) = send_child_data(j,k,2,c,h) + &
                        Py(jj,j,c)*Pz(kk,k,c)*Pew(ii,2)*work(nguard+n1off(c)+nxb/2+ii,nguard+n2off(c)+1+(j-1)/2+jj,nguard+n3off(c)+1+(k-1)/2+kk,b)
                  enddo
                enddo
              enddo
            enddo
          enddo

          ! U & D faces

          do j = 1, nyb
            do i = 1, nxb
              do jj = -2, 2
                do ii = -2, 2
                  do kk = -2, 2
                    send_child_data(i,j,5,c,h) = send_child_data(i,j,5,c,h) + &
                        Px(ii,i,c)*Py(jj,j,c)*Pud(kk,1)*work(nguard+n1off(c)+1+(i-1)/2+ii,nguard+n2off(c)+1+(j-1)/2+jj,nguard+n3off(c)+1+kk,b)
                    send_child_data(i,j,6,c,h) = send_child_data(i,j,6,c,h) + &
                        Px(ii,i,c)*Py(jj,j,c)*Pud(kk,2)*work(nguard+n1off(c)+1+(i-1)/2+ii,nguard+n2off(c)+1+(j-1)/2+jj,nguard+n3off(c)+nzb/2+kk,b)
                  enddo
                enddo
              enddo
            enddo
          enddo

        endif
      enddo

    endif

    call dBaseReleaseDataPtrSingleBlock(b, solnData)

! Next, loop over children.  If a child is on this processor, set its
! boundary values directly.  If it is off-processor, send the data to
! the owning processor (non-blocking).

    any_sent = .false.

    do c = 1, nchild
      if ((ichild == 0) .or. (ichild == c)) then
        if (child(2,c,b) == mype) then  ! local child
          solnData => dBaseGetDataPtrSingleBlock(child(1,c,b), GC)
          solnData(ito,nguard,nguard*k2d+1:nguard*k2d+nyb,&
                   nguard*k3d+1:nguard*k3d+nzb) = &
                                         send_child_data(1:nyb,1:nzb,1,c,h)
          solnData(ito,nguard+nxb+1,nguard*k2d+1:nguard*k2d+nyb,&
                   nguard*k3d+1:nguard*k3d+nzb) = &
                                         send_child_data(1:nyb,1:nzb,2,c,h)
          if (ndim >= 2) then
            solnData(ito,nguard+1:nguard+nxb,nguard,&
                     nguard*k3d+1:nguard*k3d+nzb) = &
                                         send_child_data(1:nxb,1:nzb,3,c,h)
            solnData(ito,nguard+1:nguard+nxb,nguard+nyb+1,&
                     nguard*k3d+1:nguard*k3d+nzb) = &
                                         send_child_data(1:nxb,1:nzb,4,c,h)
          endif
          if (ndim == 3) then
            solnData(ito,nguard+1:nguard+nxb,nguard+1:nguard+nyb,&
                     nguard)           = send_child_data(1:nxb,1:nyb,5,c,h)
            solnData(ito,nguard+1:nguard+nxb,nguard+1:nguard+nyb,&
                     nguard+nzb+1)     = send_child_data(1:nxb,1:nyb,6,c,h)
          endif
          call dBaseReleaseDataPtrSingleBlock(child(1,c,b), solnData)
        else                            ! remote child
          any_sent = .true.
          nsent = nsent + 1
          call mpi_issend(send_child_data(1,1,1,c,h), nmax1*nmax2*nfaces, &
                          MPI_DOUBLE_PRECISION, child(2,c,b), child(1,c,b), &
                          MPI_COMM_WORLD, send_req(nsent), ierr)
        endif
      endif
    enddo

    if (any_sent) h = h + 1
    if (h > nbbuf) &
      call abort_flash("Buffer space exceeded in hg_prolong_bndries")

  endif
enddo

do b = 1, lnblocks
  if ((dBaseRefinementLevel(b) == level+1) .and. (parent(2,b) /= mype)) then

! If parent is on another processor, receive the boundary data from
! the parent (blocking).

    solnData => dBaseGetDataPtrSingleBlock(b, GC)

    call mpi_recv(recv_child_data(1,1,1), nmax1*nmax2*nfaces, &
                   MPI_DOUBLE_PRECISION, parent(2,b), b, &
                   MPI_COMM_WORLD, status, ierr)

    solnData(ito,nguard,nguard*k2d+1:nguard*k2d+nyb,&
             nguard*k3d+1:nguard*k3d+nzb) = &
                                   recv_child_data(1:nyb,1:nzb,1)
    solnData(ito,nguard+nxb+1,nguard*k2d+1:nguard*k2d+nyb,&
             nguard*k3d+1:nguard*k3d+nzb) = &
                                   recv_child_data(1:nyb,1:nzb,2)
    if (ndim >= 2) then
      solnData(ito,nguard+1:nguard+nxb,nguard,&
               nguard*k3d+1:nguard*k3d+nzb) = &
                                   recv_child_data(1:nxb,1:nzb,3)
      solnData(ito,nguard+1:nguard+nxb,nguard+nyb+1,&
               nguard*k3d+1:nguard*k3d+nzb) = &
                                   recv_child_data(1:nxb,1:nzb,4)
    endif
    if (ndim == 3) then
      solnData(ito,nguard+1:nguard+nxb,nguard+1:nguard+nyb,&
               nguard)           = recv_child_data(1:nxb,1:nyb,5)
      solnData(ito,nguard+1:nguard+nxb,nguard+1:nguard+nyb,&
               nguard+nzb+1)     = recv_child_data(1:nxb,1:nyb,6)
    endif

    call dBaseReleaseDataPtrSingleBlock(b, solnData)

  endif
enddo

call mpi_waitall(nsent, send_req, MPI_STATUSES_IGNORE, ierr)

call timer_stop("hg_prolong_bndries")

!==============================================================================

return
end subroutine hg_prolong_bndries
