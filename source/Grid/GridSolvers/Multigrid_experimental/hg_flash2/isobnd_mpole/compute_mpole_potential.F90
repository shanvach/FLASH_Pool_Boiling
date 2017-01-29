!******************************************************************************

!  Routine:     compute_mpole_potential()
 
!  Description: Computes the potential field due to a density field associated
!               with a given set of multipole moments.  The moments are taken
!               from the mpole_common variable Moment().  On output, the
!               variable indexed by ipotvar contains the potential.  This
!               calculation is entirely local to each processor, as each
!               processor now has a separate copy of the moments.

!               Special version for isolated boundary image mass.


      subroutine compute_mpole_potential (ipotvar, poisfact)

!==============================================================================

      use mpole_common
      use physical_constants
      use dBase, ONLY: dBasePropertyInteger, dBaseBlockSize, dBaseBlockCoord, &
                       dBaseGetDataPtrSingleBlock, nguard, GC, &
                       dBaseNeighborBlockList, dBaseLocalBlockCount, &
                       dBaseReleaseDataPtrSingleBlock
      implicit none

      integer   :: ipotvar
      real      :: poisfact

      integer   :: i, j, k, lb, lnblocks
      real      :: potential, delx, dely, delz, dvol, x, y, z
      real      :: mpfactor

      integer   :: ii, jj, kk, Nint = 6
      real      :: xx, yy, zz, potsum
      real      :: delxx, delyy, delzz, ddvol

      integer, parameter :: MAXDIM = 3
      real               :: blksize(MAXDIM), blkcoord(MAXDIM), nbrs(2*MAXDIM)

      real, pointer      :: solnData(:,:,:,:)

      real, save    :: Nint_inv, fourpi_inv, nxb_inv, nyb_inv, nzb_inv
      logical, save :: FirstCall = .true.

!==============================================================================

!               One-time-only initializations.

      if (FirstCall) then

         fourpi_inv = 1.0/fourpi

         nxb_inv = 1.0/float(nxb)
         nyb_inv = 1.0/float(nyb)
         nzb_inv = 1.0/float(nzb)
         
         Nint_inv = 1.0/float(Nint)

         FirstCall = .false.

      endif

!               Scale the Poisson source term factor appropriately.

      mpfactor = poisfact * fourpi_inv

!               Compute potential on all locally held blocks.

      lnblocks = dBaseLocalBlockCount()

      do lb = 1, lnblocks

!               Get size information for this block.

        blksize  = dBaseBlockSize(lb)
        blkcoord = dBaseBlockCoord(lb) - 0.5*blksize

!               Compute dimensions of each zone and subzone.

        delx = blksize(1) * nxb_inv
        dely = blksize(2) * nyb_inv * k2d
        delz = blksize(3) * nzb_inv * k3d

        delxx = delx * Nint_inv
        delyy = dely * Nint_inv
        delzz = delz * Nint_inv

!               Get pointer to solution data.

        solnData => dBaseGetDataPtrSingleBlock(lb, GC)

        nbrs = dBaseNeighborBlockList(lb)

!               Compute potential boundary values.

      select case (mpole_geometry)

        case (G_3DCARTESIAN)

          if (nbrs(1) <= -20) then
            xx = blkcoord(1) - Xcm
            do k = kmin, kmax
              zz = blkcoord(3) + (k-kmin+0.5)*delz - Zcm
              do j = jmin, jmax
                yy = blkcoord(2) + (j-jmin+0.5)*dely - Ycm
                call zone_potential (xx, yy, zz, potential)
                solnData(ipotvar,nguard,j,k) = -mpfactor*potential
              enddo
            enddo
          endif

          if (nbrs(2) <= -20) then
            xx = blkcoord(1) + blksize(1) - Xcm
            do k = kmin, kmax
              zz = blkcoord(3) + (k-kmin+0.5)*delz - Zcm
              do j = jmin, jmax
                yy = blkcoord(2) + (j-jmin+0.5)*dely - Ycm
                call zone_potential (xx, yy, zz, potential)
                solnData(ipotvar,nguard+nxb+1,j,k) = -mpfactor*potential
              enddo
            enddo
          endif

          if (nbrs(3) <= -20) then
            yy = blkcoord(2) - Ycm
            do k = kmin, kmax
              zz = blkcoord(3) + (k-kmin+0.5)*delz - Zcm
              do i = imin, imax
                xx = blkcoord(1) + (i-imin+0.5)*delx - Xcm
                call zone_potential (xx, yy, zz, potential)
                solnData(ipotvar,i,nguard,k) = -mpfactor*potential
              enddo
            enddo
          endif

          if (nbrs(4) <= -20) then
            yy = blkcoord(2) + blksize(2) - Ycm
            do k = kmin, kmax
              zz = blkcoord(3) + (k-kmin+0.5)*delz - Zcm
              do i = imin, imax
                xx = blkcoord(1) + (i-imin+0.5)*delx - Xcm
                call zone_potential (xx, yy, zz, potential)
                solnData(ipotvar,i,nguard+nyb+1,k) = -mpfactor*potential
              enddo
            enddo
          endif

          if (nbrs(5) <= -20) then
            zz = blkcoord(3) - Zcm
            do j = jmin, jmax
              yy = blkcoord(2) + (j-jmin+0.5)*dely - Ycm
              do i = imin, imax
                xx = blkcoord(1) + (i-imin+0.5)*delx - Xcm
                call zone_potential (xx, yy, zz, potential)
                solnData(ipotvar,i,j,nguard) = -mpfactor*potential
              enddo
            enddo
          endif

          if (nbrs(6) <= -20) then
            zz = blkcoord(3) + blksize(3) - Zcm
            do j = jmin, jmax
              yy = blkcoord(2) + (j-jmin+0.5)*dely - Ycm
              do i = imin, imax
                xx = blkcoord(1) + (i-imin+0.5)*delx - Xcm
                call zone_potential (xx, yy, zz, potential)
                solnData(ipotvar,i,j,nguard+nzb+1) = -mpfactor*potential
              enddo
            enddo
          endif

        case (G_2DCYLINDRICAL)

          if (nbrs(2) <= -20) then
            xx = blkcoord(1) + blksize(1) - Xcm
            zz = 0.
            do j = jmin, jmax
              yy = blkcoord(2) + (j-jmin+0.5)*dely - Ycm
              call zone_potential (xx, 0., yy, potential)
              solnData(ipotvar,nguard+nxb+1,j,1) = -mpfactor*potential
            enddo
          endif

          if (nbrs(3) <= -20) then
            yy = blkcoord(2) - Ycm
            zz = 0.
            do i = imin, imax
              xx = blkcoord(1) + (i-imin+0.5)*delx - Xcm
              call zone_potential (xx, 0., yy, potential)
              solnData(ipotvar,i,nguard,1) = -mpfactor*potential
            enddo
          endif

          if (nbrs(4) <= -20) then
            yy = blkcoord(2) + blksize(2) - Ycm
            zz = 0.
            do i = imin, imax
              xx = blkcoord(1) + (i-imin+0.5)*delx - Xcm
              call zone_potential (xx, 0., yy, potential)
              solnData(ipotvar,i,nguard+nyb+1,1) = -mpfactor*potential
            enddo
          endif

        case (G_1DSPHERICAL)

          if (nbrs(2) <= -20) then
            xx = blkcoord(1) + blksize(1) - Xcm
            yy = 0.
            zz = 0.
            call zone_potential (xx, 0., 0., potential)
            solnData(ipotvar,nguard+nxb+1,1,1) = -mpfactor*potential
          endif

      end select

        call dBaseReleaseDataPtrSingleBlock (lb, solnData)

      enddo

!==============================================================================

      return
      end
