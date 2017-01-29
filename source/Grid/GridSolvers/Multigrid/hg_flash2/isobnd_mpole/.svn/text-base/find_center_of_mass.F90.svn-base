!******************************************************************************

!  Routine:     find_center_of_mass()
 
!  Description: Computes the centroid of the specified variable and returns
!               its location in the mpole_common variables Xcm, Ycm, and Zcm.
!               Also computes the total value of the quantity and leaves it
!               in the variable Mtot.  If Mtot=0, the routine aborts.

!               Special version for isolated boundary image mass.


      subroutine find_center_of_mass (idensvar)

!==============================================================================

      use mpole_common
      use dBase, ONLY: dBaseNodeType, dBaseLocalBlockCount
      use logfile
      implicit none

      include "mpif.h"

      integer :: idensvar

      integer, parameter :: nsum = 4    ! Number of summed quantities;
                                        ! coordinate w/mpole_sum_local
      real    :: sum(nsum), lsum(nsum)

      integer :: i, error, lnblocks

!==============================================================================

!               Sum quantities over all locally held leaf blocks.

      sum  = 0.
      lsum = 0.

      lnblocks = dBaseLocalBlockCount()

      do i = 1, lnblocks
        if (dBaseNodeType(i) == 1) &
          call mpole_sum_local (lsum, nsum, i, idensvar)
      enddo

!==============================================================================

!               Give all processors a copy of the global sums.

      call mpi_allreduce (lsum, sum, nsum, MPI_DOUBLE_PRECISION, & 
                          MPI_SUM, MPI_COMM_WORLD, error)

!               Now normalize the center-of-mass coordinates.

      Mtot = sum(1)
      if (abs(Mtot) < 1.E-99) &
        call abort_flash ("FATAL:  find_center_of_mass:  Mtot = 0")
      Xcm = sum(2) / Mtot
      Ycm = sum(3) / Mtot
      Zcm = sum(4) / Mtot

!==============================================================================

      return
      end


!******************************************************************************

!  Routine:     mpole_sum_local()

!  Description: Accumulate values of the density and density-weighted position
!               from the interior of a given block.  The results are stored in
!               a given vector (lsum), which is assumed to have been
!               initialized by the calling routine.


      subroutine mpole_sum_local (lsum, nsum, block_no, idensvar)

!==============================================================================

      use mpole_common
      use dBase, ONLY: dBaseBlockSize, dBaseBlockCoord, &
                       dBaseGetDataPtrSingleBlock,      &
                       dBaseReleaseDataPtrSingleBlock,  &
                       dBaseNeighborBlockList, nguard, GC
      implicit none

      integer :: block_no, idensvar, nsum
      real    :: lsum(nsum)

      integer, parameter :: MAXDIM = 3
      real               :: size(MAXDIM), coord(MAXDIM)
      integer            :: nbrs(2*MAXDIM)
      real               :: xx, yy, zz
      real, pointer      :: solnData(:,:,:,:)
      integer            :: i, j, k, n
      real               :: dvol, delx, dely, delz, delm

      real               :: cylfactor

      logical, save      :: firstCall = .true.

!==============================================================================

! if we are doing only a quadrant of the star in 2-d cylindrical coords, then
! we pick up an extra factor of 2 in the volume calculations because of the 
! symmetry about the y = 0 axis

      if (firstCall) then

        if (quadrant) then
          cylfactor = 2.0*twopi
        else
          cylfactor = twopi
        endif

        firstCall = .false.

      endif

!==============================================================================

!               Get size information for this block.

      size  = dBaseBlockSize(block_no)
      coord = dBaseBlockCoord(block_no) - 0.5*size

!               Compute dimensions of each zone.

      delx = size(1) / real(nxb)
      if (ndim >= 2) dely = size(2) / real(nyb)
      if (ndim == 3) delz = size(3) / real(nzb)

!               Obtain a pointer to the block and its neighbor list.

      solnData => dBaseGetDataPtrSingleBlock(block_no, GC)
      nbrs = dBaseNeighborBlockList(block_no)

!               Sum contributions from this block.


      select case (mpole_geometry)

        case (G_3DCARTESIAN)

          if (nbrs(1) <= -20) then
            xx = coord(1)
            do k = kmin, kmax
              zz = coord(3) + (k-kmin+0.5)*delz
              do j = jmin, jmax
                yy = coord(2) + (j-jmin+0.5)*dely
                delm    = solnData(idensvar,nguard,j,k) * dely * delz
                lsum(1) = lsum(1) + delm
                lsum(2) = lsum(2) + delm*xx
                lsum(3) = lsum(3) + delm*yy
                lsum(4) = lsum(4) + delm*zz
              enddo
            enddo
          endif

          if (nbrs(2) <= -20) then
            xx = coord(1) + size(1)
            do k = kmin, kmax
              zz = coord(3) + (k-kmin+0.5)*delz
              do j = jmin, jmax
                yy = coord(2) + (j-jmin+0.5)*dely
                delm    = solnData(idensvar,nguard+nxb+1,j,k) * dely * delz
                lsum(1) = lsum(1) + delm
                lsum(2) = lsum(2) + delm*xx
                lsum(3) = lsum(3) + delm*yy
                lsum(4) = lsum(4) + delm*zz
              enddo
            enddo
          endif

          if (nbrs(3) <= -20) then
            yy = coord(2)
            do k = kmin, kmax
              zz = coord(3) + (k-kmin+0.5)*delz
              do i = imin, imax
                xx = coord(1) + (i-imin+0.5)*delx
                delm    = solnData(idensvar,i,nguard,k) * delx * delz
                lsum(1) = lsum(1) + delm
                lsum(2) = lsum(2) + delm*xx
                lsum(3) = lsum(3) + delm*yy
                lsum(4) = lsum(4) + delm*zz
              enddo
            enddo
          endif

          if (nbrs(4) <= -20) then
            yy = coord(2) + size(2)
            do k = kmin, kmax
              zz = coord(3) + (k-kmin+0.5)*delz
              do i = imin, imax
                xx = coord(1) + (i-imin+0.5)*delx
                delm    = solnData(idensvar,i,nguard+nyb+1,k) * delx * delz
                lsum(1) = lsum(1) + delm
                lsum(2) = lsum(2) + delm*xx
                lsum(3) = lsum(3) + delm*yy
                lsum(4) = lsum(4) + delm*zz
              enddo
            enddo
          endif

          if (nbrs(5) <= -20) then
            zz = coord(3)
            do j = jmin, jmax
              yy = coord(2) + (j-jmin+0.5)*dely
              do i = imin, imax
                xx = coord(1) + (i-imin+0.5)*delx
                delm    = solnData(idensvar,i,j,nguard) * delx * dely
                lsum(1) = lsum(1) + delm
                lsum(2) = lsum(2) + delm*xx
                lsum(3) = lsum(3) + delm*yy
                lsum(4) = lsum(4) + delm*zz
              enddo
            enddo
          endif

          if (nbrs(6) <= -20) then
            zz = coord(3) + size(3)
            do j = jmin, jmax
              yy = coord(2) + (j-jmin+0.5)*dely
              do i = imin, imax
                xx = coord(1) + (i-imin+0.5)*delx
                delm    = solnData(idensvar,i,j,nguard+nzb+1) * delx * dely
                lsum(1) = lsum(1) + delm
                lsum(2) = lsum(2) + delm*xx
                lsum(3) = lsum(3) + delm*yy
                lsum(4) = lsum(4) + delm*zz
              enddo
            enddo
          endif

        case (G_2DCYLINDRICAL)

          if (nbrs(2) <= -20) then
            xx = coord(1) + size(1)
            do j = jmin, jmax
              yy = coord(2) + (j-jmin+0.5)*dely
              delm = solnData(idensvar,nguard+nxb+1,j,1) * cylfactor * xx * dely
              lsum(1) = lsum(1) + delm
              ! center of mass is on x-axis in 2D cylindrical
              ! center of mass is on y-axis if doing only a quadrant
              if (.not. quadrant) lsum(3) = lsum(3) + delm*yy
            enddo
          endif

          if (nbrs(3) <= -20) then
            yy = coord(2)
            do i = imin, imax
              xx = coord(1) + (i-imin+0.5)*delx
              delm = solnData(idensvar,i,nguard,1) * cylfactor * xx * delx
              lsum(1) = lsum(1) + delm
              ! center of mass is on x-axis in 2D cylindrical
              ! center of mass is on y-axis if doing only a quadrant
              if (.not. quadrant) lsum(3) = lsum(3) + delm*yy
            enddo
          endif

          if (nbrs(4) <= -20) then
            yy = coord(2) + size(2)
            do i = imin, imax
              xx = coord(1) + (i-imin+0.5)*delx
              delm = solnData(idensvar,i,nguard,1) * cylfactor * xx * delx
              lsum(1) = lsum(1) + delm
              ! center of mass is on x-axis in 2D cylindrical
              ! center of mass is on y-axis if doing only a quadrant
              if (.not. quadrant) lsum(3) = lsum(3) + delm*yy
            enddo
          endif

        case (G_1DSPHERICAL)

          if (nbrs(2) <= -20) then
            xx = coord(1) + size(1)
            delm = solnData(idensvar,nguard+nxb+1,1,1) * fourpi * xx**2
            lsum(1) = lsum(1) + delm
            ! center of mass is at origin in 1D spherical
          endif

      end select

      call dBaseReleaseDataPtrSingleBlock (block_no, solnData)

!==============================================================================

      return
      end
