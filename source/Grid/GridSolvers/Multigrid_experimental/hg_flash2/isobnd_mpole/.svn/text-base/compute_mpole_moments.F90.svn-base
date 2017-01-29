!!****if* source/Grid/GridSolvers/Multigrid_experimental/hg_flash2/isobnd_mpole/compute_mpole_moments
!!
!! NAME
!!
!!  compute_mpole_moments
!!
!! 
!! SYNOPSIS
!!
!!  compute_mpole_moments(idensvar)
!!
!!  compute_mpole_moments(integer)
!!
!!
!! DESCRIPTION
!!
!!  Compute the multipole moments of the density distribution,
!!  assuming the center of mass and the total mass have first been
!!  computed by find_center_of_mass().  On output, the Moment()
!!  array declared in mpole_common.F90 contains the moments.
!!
!!  Special version for isolated boundary image mass.
!!
!!***

      subroutine compute_mpole_moments (idensvar)

!==============================================================================

      use mpole_common
      use dBase, ONLY: dBaseNodeType, dBasePropertyInteger, &
                       dBaseGetDataPtrSingleBlock, dBaseBlockSize, &
                       dBaseBlockCoord, dBaseNeighborBlockList, nguard, GC, &
                       dBaseLocalBlockCount, dBaseReleaseDataPtrSingleBlock
      implicit none
      include "mpif.h"


      integer   :: idensvar

      integer   :: q, i, j, k, l, m, error, lnblocks
      real      :: zonemass, dvol, delx, dely, delz, x, y, z

      integer   :: lb, ii, jj, kk, Nint = 2
      real      :: xx, yy, zz
      real      :: delxx, delyy, delzz
      real      :: zonedens

      integer, parameter :: MAXDIM = 3
      real               :: size(MAXDIM), coord(MAXDIM)
      integer            :: nbrs(2*MAXDIM)

      real, pointer      :: solnData(:,:,:,:)

      real, save         :: cylfactor

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

!               Sum quantities over all locally held leaf blocks.

      Moment(:,:,:,:,:) = 0.

      lnblocks = dBaseLocalBlockCount()

      do lb = 1, lnblocks

        if (dBaseNodeType(lb) == 1) then        ! Is a leaf block.

!               Get size information for this block.

          size  = dBaseBlockSize(lb)
          coord = dBaseBlockCoord(lb) - 0.5*size

!               Compute dimensions of each zone and subzone.

          delx = size(1) / real(nxb)
          dely = size(2) / real(nyb) * k2d
          delz = size(3) / real(nzb) * k3d

          delxx = delx / real(Nint)
          delyy = dely / real(Nint)
          delzz = delz / real(Nint)

!               Get pointer to solution data.

          solnData => dBaseGetDataPtrSingleBlock(lb, GC)

      nbrs = dBaseNeighborBlockList(lb)

!               Compute the moment contributions for this block.

      select case (mpole_geometry)

        case (G_3DCARTESIAN)

          if (nbrs(1) <= -20) then
            xx = coord(1) - Xcm
            do k = kmin, kmax
              zz = coord(3) + (k-kmin+0.5)*delz - Zcm
              do j = jmin, jmax
                yy = coord(2) + (j-jmin+0.5)*dely - Ycm
                zonemass = solnData(idensvar,nguard,j,k) * dely * delz
                call zone_moments (xx, yy, zz, zonemass)
              enddo
            enddo
          endif

          if (nbrs(2) <= -20) then
            xx = coord(1) + size(1) - Xcm
            do k = kmin, kmax
              zz = coord(3) + (k-kmin+0.5)*delz - Zcm
              do j = jmin, jmax
                yy = coord(2) + (j-jmin+0.5)*dely - Ycm
                zonemass = solnData(idensvar,nguard+nxb+1,j,k) * dely * delz
                call zone_moments (xx, yy, zz, zonemass)
              enddo
            enddo
          endif

          if (nbrs(3) <= -20) then
            yy = coord(2) - Ycm
            do k = kmin, kmax
              zz = coord(3) + (k-kmin+0.5)*delz - Zcm
              do i = imin, imax
                xx = coord(1) + (i-imin+0.5)*delx - Xcm
                zonemass = solnData(idensvar,i,nguard,k) * delx * delz
                call zone_moments (xx, yy, zz, zonemass)
              enddo
            enddo
          endif

          if (nbrs(4) <= -20) then
            yy = coord(2) + size(2) - Ycm
            do k = kmin, kmax
              zz = coord(3) + (k-kmin+0.5)*delz - Zcm
              do i = imin, imax
                xx = coord(1) + (i-imin+0.5)*delx - Xcm
                zonemass = solnData(idensvar,i,nguard+nyb+1,k) * delx * delz
                call zone_moments (xx, yy, zz, zonemass)
              enddo
            enddo
          endif

          if (nbrs(5) <= -20) then
            zz = coord(3) - Zcm
            do j = jmin, jmax
              yy = coord(2) + (j-jmin+0.5)*dely - Ycm
              do i = imin, imax
                xx = coord(1) + (i-imin+0.5)*delx - Xcm
                zonemass = solnData(idensvar,i,j,nguard) * delx * dely
                call zone_moments (xx, yy, zz, zonemass)
              enddo
            enddo
          endif

          if (nbrs(6) <= -20) then
            zz = coord(3) + size(3) - Zcm
            do j = jmin, jmax
              yy = coord(2) + (j-jmin+0.5)*dely - Ycm
              do i = imin, imax
                xx = coord(1) + (i-imin+0.5)*delx - Xcm
                zonemass = solnData(idensvar,i,j,nguard+nzb+1) * delx * dely
                call zone_moments (xx, yy, zz, zonemass)
              enddo
            enddo
          endif

        case (G_2DCYLINDRICAL)

          if (nbrs(2) <= -20) then
            xx = coord(1) + size(1) - Xcm
            zz = 0.
            do j = jmin, jmax
              yy = coord(2) + (j-jmin+0.5)*dely - Ycm
              zonemass = solnData(idensvar,nguard+nxb+1,j,1) * &
                         cylfactor * xx * dely
              call zone_moments (xx, 0., yy, zonemass)
            enddo
          endif

          if (nbrs(3) <= -20) then
            yy = coord(2) - Ycm
            zz = 0.
            do i = imin, imax
              xx = coord(1) + (i-imin+0.5)*delx - Xcm
              zonemass = solnData(idensvar,i,nguard,1) * &
                         cylfactor * xx * delx
              call zone_moments (xx, 0., yy, zonemass)
            enddo
          endif

          if (nbrs(4) <= -20) then
            yy = coord(2) + size(2) - Ycm
            zz = 0.
            do i = imin, imax
              xx = coord(1) + (i-imin+0.5)*delx - Xcm
              zonemass = solnData(idensvar,i,nguard+nyb+1,1) * &
                         cylfactor * xx * delx
              call zone_moments (xx, 0., yy, zonemass)
            enddo
          endif

        case (G_1DSPHERICAL)

          if (nbrs(2) <= -20) then
            xx = coord(1) + size(1) - Xcm
            yy = 0.
            zz = 0.
            zonemass = solnData(idensvar,nguard+nxb+1,1,1) * fourpi * xx**2
            call zone_moments (xx, 0., 0., zonemass)
          endif

       end select

       call dBaseReleaseDataPtrSingleBlock (lb, solnData)

       endif

     enddo

! In zone_moments, we only added the contribution of each zone to its
! particular radial bin.  Each zone should contribute its inner moment to
! all zones with radius greater than it, and its outer moment to all zones
! with radius less than it.

     do m = 0, mpole_lmax
       do l = m, mpole_lmax

             do i = Even, Odd

                do q = 2, qmax
                   Moment(q,i,Inner,l,m) = Moment(q,i,Inner,l,m) + &
                        Moment(q-1,i,Inner,l,m)
                enddo

                do q = qmax-1, 1, -1
                   Moment(q,i,Outer,l,m) = Moment(q,i,Outer,l,m) + &
                        Moment(q+1,i,Outer,l,m)
                enddo

             enddo

        enddo
      enddo

!==============================================================================

!               Give all processors a copy of the global sums.

      do m = 0, mpole_lmax
        do l = m, mpole_lmax
          do j = Inner, Outer
            do i = Even, Odd
              call MPI_AllReduce (Moment(0,i,j,l,m), Momtmp(0), qmax+1, & 
     &                            MPI_Double_Precision, MPI_Sum,  & 
     &                            MPI_Comm_World, error)
              Moment(:,i,j,l,m) = Momtmp(:)
            enddo
          enddo
        enddo
      enddo

!               Normalize the moments properly.

      do l = 1, mpole_lmax
        do m = 1, l
          Moment(:,:,:,l,m) = Moment(:,:,:,l,m) * Leg_fact(l,m)
        enddo
      enddo

!==============================================================================

      return
      end
