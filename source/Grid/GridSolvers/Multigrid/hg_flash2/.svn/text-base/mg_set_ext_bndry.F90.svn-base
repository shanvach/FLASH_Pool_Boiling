!******************************************************************************* 

!  Routine:     mg_set_ext_bndry()
 
!  Description: Set values of fluid variables in the ghost zones of all blocks.
!               The type of boundary for a given block in a given direction is
!               specified by the neigh() array:  if positive, this gives the
!               block number of the neighboring block (and the number of the
!               processor owning that block); if negative, the boundary is
!               external, and the boundary value indicates the type of boundary
!               to implement.  See definitions.fh for the definitions of the
!               constants which indicate the various types of boundary.

!               Note that periodic boundaries are taken care of in the
!               initialization of neigh() and thus do not require a special
!               case here.  Also, isolated boundaries are taken care of by
!               the routine which calls the multigrid solver (e.g. poisson())
!               and thus do not require a special case here.

!               Note that for cylindrical or spherical coordinates, r=0 is
!               always a Neumann boundary.  Also, if only a single quadrant
!               is being simulated in 2D axisymmetric cylindrical coordinates,
!               z=0 is also a Neumann boundary.

!               This routine differs from tot_bnd() in that it sets ghost zones
!               for the work() array rather than the unk() array.

!  Parameters:  idiag           Whether to update corner boundaries -- seems
!                                 to be ignored!
!               idir            Direction to update:  either sweep_all (for
!                                 all directions) or sweep_#, where #=x, y, z
!               extrap          If nonzero, use boundary extrapolation rather
!                                than homogeneous Dirichlet boundaries for the
!                                exterior of the mesh.

!  Modifies:    On output, this routine should leave the ghost zones of the
!               work() array updated.  The ghost zones are zones
!               1..nguard_work for the - direction of each coordinate and
!               n#b+nguard_work+1..n#b+2*nguard_work for the + direction of
!               each coordinate, where #=x, y, z.


subroutine mg_set_ext_bndry (idiag, idir, extrap)

!===============================================================================

use mg_common
use dBase, ONLY: dBaseNeighborBlockList, dBasePropertyInteger, ndim
use dBaseDeclarations, ONLY: work, lrefine
implicit none

integer :: idiag, idir, extrap

integer :: lb, lnblocks
integer :: nbr_blks(6)

!===============================================================================

lnblocks = dBasePropertyInteger("LocalNumberOfBlocks")

! Directly set any external Dirichlet or Neumann boundaries.

if ((mg_bnd_cond == MG_BND_DIRICHLET) .or. &
    (mg_bnd_cond == MG_BND_GIVENVAL)) then              ! Dirichlet/given-value

  do lb = 1, lnblocks

    nbr_blks = dBaseNeighborBlockList(lb)

    if ((mg_geometry == MG_GEOM_1DCARTESIAN) .or. &
        (mg_geometry == MG_GEOM_2DCARTESIAN) .or. &
        (mg_geometry == MG_GEOM_3DCARTESIAN)) then

      if (nbr_blks(1) <= -20) then
        if (extrap /= 0) then
          work(ili-1,:,:,lb) = 2.*work(ili,:,:,lb) - work(ili+1,:,:,lb)
          work(ili-2,:,:,lb) = 3.*work(ili,:,:,lb) - 2.*work(ili+1,:,:,lb)
          work(ili-3,:,:,lb) = 4.*work(ili,:,:,lb) - 3.*work(ili+1,:,:,lb)
          work(ili-4,:,:,lb) = 5.*work(ili,:,:,lb) - 4.*work(ili+1,:,:,lb)
        else
          work(ili-1,:,:,lb) = -   work(ili,:,:,lb)
          work(ili-2,:,:,lb) = -3.*work(ili,:,:,lb)
          work(ili-3,:,:,lb) = -5.*work(ili,:,:,lb)
          work(ili-4,:,:,lb) = -7.*work(ili,:,:,lb)
        endif
      endif

    else

      if (nbr_blks(1) <= -20) work(ili-1,:,:,lb) = work(ili,:,:,lb)

    endif

    if (nbr_blks(2) <= -20) then
      if (extrap /= 0) then
        work(iui+1,:,:,lb) = 2.*work(iui,:,:,lb) - work(iui-1,:,:,lb)
        work(iui+2,:,:,lb) = 3.*work(iui,:,:,lb) - 2.*work(iui-1,:,:,lb)
        work(iui+3,:,:,lb) = 4.*work(iui,:,:,lb) - 3.*work(iui-1,:,:,lb)
        work(iui+4,:,:,lb) = 5.*work(iui,:,:,lb) - 4.*work(iui-1,:,:,lb)
      else
        work(iui+1,:,:,lb) = -work(iui,:,:,lb)
        work(iui+2,:,:,lb) = -3.*work(iui,:,:,lb)
        work(iui+3,:,:,lb) = -5.*work(iui,:,:,lb)
        work(iui+4,:,:,lb) = -7.*work(iui,:,:,lb)
      endif
    endif

    if (ndim >= 2) then

      if ((mg_geometry /= MG_GEOM_2DCYLAXISYM) .or. (.not. quadrant)) then
        if (nbr_blks(3) <= -20) then
          if (extrap /= 0) then
            work(:,jli-1,:,lb) = 2.*work(:,jli,:,lb) - work(:,jli+1,:,lb)
            work(:,jli-2,:,lb) = 3.*work(:,jli,:,lb) - 2.*work(:,jli+1,:,lb)
            work(:,jli-3,:,lb) = 4.*work(:,jli,:,lb) - 3.*work(:,jli+1,:,lb)
            work(:,jli-4,:,lb) = 5.*work(:,jli,:,lb) - 4.*work(:,jli+1,:,lb)
          else
            work(:,jli-1,:,lb) = -work(:,jli,:,lb)
            work(:,jli-2,:,lb) = -3.*work(:,jli,:,lb)
            work(:,jli-3,:,lb) = -5.*work(:,jli,:,lb)
            work(:,jli-4,:,lb) = -7.*work(:,jli,:,lb)
          endif
        endif
        if ((nbr_blks(1) <= -20) .and. (nbr_blks(3) <= -20)) then
          work(ili-1,jli-1,:,lb) =  10.*work(ili,jli,:,lb)     - 10.*work(ili,jli+1,:,lb)   +  5.*work(ili,jli+2,:,lb)   -      work(ili,jli+3,:,lb) &
                                  - 10.*work(ili+1,jli,:,lb)   +  5.*work(ili+1,jli+1,:,lb) -     work(ili+1,jli+2,:,lb) +  5.* work(ili+2,jli,:,lb) &
                                  -     work(ili+2,jli+1,:,lb) -     work(ili+3,jli,:,lb)
          work(ili-1,jli-2,:,lb) =  20.*work(ili,jli,:,lb)     - 30.*work(ili,jli+1,:,lb)   + 18.*work(ili,jli+2,:,lb)   -  4.* work(ili,jli+3,:,lb) &
                                  - 15.*work(ili+1,jli,:,lb)   + 12.*work(ili+1,jli+1,:,lb) -  3.*work(ili+1,jli+2,:,lb) +  6.* work(ili+2,jli,:,lb) &
                                  -  2.*work(ili+2,jli+1,:,lb) -     work(ili+3,jli,:,lb)
          work(ili-2,jli-1,:,lb) =  20.*work(ili,jli,:,lb)     - 15.*work(ili,jli+1,:,lb)   +  6.*work(ili,jli+2,:,lb)   -      work(ili,jli+3,:,lb) &
                                  - 30.*work(ili+1,jli,:,lb)   + 12.*work(ili+1,jli+1,:,lb) -  2.*work(ili+1,jli+2,:,lb) + 18.* work(ili+2,jli,:,lb) &
                                  -  3.*work(ili+2,jli+1,:,lb) -  4.*work(ili+3,jli,:,lb)
          work(ili-2,jli-2,:,lb) =  35.*work(ili,jli,:,lb)     - 42.*work(ili,jli+1,:,lb)   + 21.*work(ili,jli+2,:,lb)   -  4.* work(ili,jli+3,:,lb) &
                                  - 42.*work(ili+1,jli,:,lb)   + 28.*work(ili+1,jli+1,:,lb) -  6.*work(ili+1,jli+2,:,lb) + 21.* work(ili+2,jli,:,lb) &
                                  -  6.*work(ili+2,jli+1,:,lb) -  4.*work(ili+3,jli,:,lb)
        endif
        if ((nbr_blks(2) <= -20) .and. (nbr_blks(3) <= -20)) then
          work(iui+1,jli-1,:,lb) =  10.*work(iui,jli,:,lb)     - 10.*work(iui,jli+1,:,lb)   +  5.*work(iui,jli+2,:,lb)   -      work(iui,jli+3,:,lb) &
                                  - 10.*work(iui-1,jli,:,lb)   +  5.*work(iui-1,jli+1,:,lb) -     work(iui-1,jli+2,:,lb) +  5.* work(iui-2,jli,:,lb) &
                                  -     work(iui-2,jli+1,:,lb) -     work(iui-3,jli,:,lb)
          work(iui+1,jli-2,:,lb) =  20.*work(iui,jli,:,lb)     - 30.*work(iui,jli+1,:,lb)   + 18.*work(iui,jli+2,:,lb)   -  4.* work(iui,jli+3,:,lb) &
                                  - 15.*work(iui-1,jli,:,lb)   + 12.*work(iui-1,jli+1,:,lb) -  3.*work(iui-1,jli+2,:,lb) +  6.* work(iui-2,jli,:,lb) &
                                  -  2.*work(iui-2,jli+1,:,lb) -     work(iui-3,jli,:,lb)
          work(iui+2,jli-1,:,lb) =  20.*work(iui,jli,:,lb)     - 15.*work(iui,jli+1,:,lb)   +  6.*work(iui,jli+2,:,lb)   -      work(iui,jli+3,:,lb) &
                                  - 30.*work(iui-1,jli,:,lb)   + 12.*work(iui-1,jli+1,:,lb) -  2.*work(iui-1,jli+2,:,lb) + 18.* work(iui-2,jli,:,lb) &
                                  -  3.*work(iui-2,jli+1,:,lb) -  4.*work(iui-3,jli,:,lb)
          work(iui+2,jli-2,:,lb) =  35.*work(iui,jli,:,lb)     - 42.*work(iui,jli+1,:,lb)   + 21.*work(iui,jli+2,:,lb)   -  4.* work(iui,jli+3,:,lb) &
                                  - 42.*work(iui-1,jli,:,lb)   + 28.*work(iui-1,jli+1,:,lb) -  6.*work(iui-1,jli+2,:,lb) + 21.* work(iui-2,jli,:,lb) &
                                  -  6.*work(iui-2,jli+1,:,lb) -  4.*work(iui-3,jli,:,lb)
        endif
      else
        if (nbr_blks(3) <= -20) work(:,jli-1,:,lb) = work(:,jli,:,lb)
      endif

      if (nbr_blks(4) <= -20) then
        if (extrap /= 0) then
          work(:,jui+1,:,lb) = 2.*work(:,jui,:,lb) - work(:,jui-1,:,lb)
          work(:,jui+2,:,lb) = 3.*work(:,jui,:,lb) - 2.*work(:,jui-1,:,lb)
          work(:,jui+3,:,lb) = 4.*work(:,jui,:,lb) - 3.*work(:,jui-1,:,lb)
          work(:,jui+4,:,lb) = 5.*work(:,jui,:,lb) - 4.*work(:,jui-1,:,lb)
        else
          work(:,jui+1,:,lb) = -work(:,jui,:,lb)
          work(:,jui+2,:,lb) = -3.*work(:,jui,:,lb)
          work(:,jui+3,:,lb) = -5.*work(:,jui,:,lb)
          work(:,jui+4,:,lb) = -7.*work(:,jui,:,lb)
        endif
      endif
      if ((nbr_blks(1) <= -20) .and. (nbr_blks(4) <= -20)) then
        work(ili-1,jui+1,:,lb) =  10.*work(ili,jui,:,lb)     - 10.*work(ili,jui-1,:,lb)   +  5.*work(ili,jui-2,:,lb)   -      work(ili,jui-3,:,lb) &
                                - 10.*work(ili+1,jui,:,lb)   +  5.*work(ili+1,jui-1,:,lb) -     work(ili+1,jui-2,:,lb) +  5.* work(ili+2,jui,:,lb) &
                                -     work(ili+2,jui-1,:,lb) -     work(ili+3,jui,:,lb)
        work(ili-1,jui+2,:,lb) =  20.*work(ili,jui,:,lb)     - 30.*work(ili,jui-1,:,lb)   + 18.*work(ili,jui-2,:,lb)   -  4.* work(ili,jui-3,:,lb) &
                                - 15.*work(ili+1,jui,:,lb)   + 12.*work(ili+1,jui-1,:,lb) -  3.*work(ili+1,jui-2,:,lb) +  6.* work(ili+2,jui,:,lb) &
                                -  2.*work(ili+2,jui-1,:,lb) -     work(ili+3,jui,:,lb)
        work(ili-2,jui+1,:,lb) =  20.*work(ili,jui,:,lb)     - 15.*work(ili,jui-1,:,lb)   +  6.*work(ili,jui-2,:,lb)   -      work(ili,jui-3,:,lb) &
                                - 30.*work(ili+1,jui,:,lb)   + 12.*work(ili+1,jui-1,:,lb) -  2.*work(ili+1,jui-2,:,lb) + 18.* work(ili+2,jui,:,lb) &
                                -  3.*work(ili+2,jui-1,:,lb) -  4.*work(ili+3,jui,:,lb)
        work(ili-2,jui+2,:,lb) =  35.*work(ili,jui,:,lb)     - 42.*work(ili,jui-1,:,lb)   + 21.*work(ili,jui-2,:,lb)   -  4.* work(ili,jui-3,:,lb) &
                                - 42.*work(ili+1,jui,:,lb)   + 28.*work(ili+1,jui-1,:,lb) -  6.*work(ili+1,jui-2,:,lb) + 21.* work(ili+2,jui,:,lb) &
                                -  6.*work(ili+2,jui-1,:,lb) -  4.*work(ili+3,jui,:,lb)
      endif
      if ((nbr_blks(2) <= -20) .and. (nbr_blks(4) <= -20)) then
        work(iui+1,jui+1,:,lb) =  10.*work(iui,jui,:,lb)     - 10.*work(iui,jui-1,:,lb)   +  5.*work(iui,jui-2,:,lb)   -      work(iui,jui-3,:,lb) &
                                - 10.*work(iui-1,jui,:,lb)   +  5.*work(iui-1,jui-1,:,lb) -     work(iui-1,jui-2,:,lb) +  5.* work(iui-2,jui,:,lb) &
                                -     work(iui-2,jui-1,:,lb) -     work(iui-3,jui,:,lb)
        work(iui+1,jui+2,:,lb) =  20.*work(iui,jui,:,lb)     - 30.*work(iui,jui-1,:,lb)   + 18.*work(iui,jui-2,:,lb)   -  4.* work(iui,jui-3,:,lb) &
                                - 15.*work(iui-1,jui,:,lb)   + 12.*work(iui-1,jui-1,:,lb) -  3.*work(iui-1,jui-2,:,lb) +  6.* work(iui-2,jui,:,lb) &
                                -  2.*work(iui-2,jui-1,:,lb) -     work(iui-3,jui,:,lb)
        work(iui+2,jui+1,:,lb) =  20.*work(iui,jui,:,lb)     - 15.*work(iui,jui-1,:,lb)   +  6.*work(iui,jui-2,:,lb)   -      work(iui,jui-3,:,lb) &
                                - 30.*work(iui-1,jui,:,lb)   + 12.*work(iui-1,jui-1,:,lb) -  2.*work(iui-1,jui-2,:,lb) + 18.* work(iui-2,jui,:,lb) &
                                -  3.*work(iui-2,jui-1,:,lb) -  4.*work(iui-3,jui,:,lb)
        work(iui+2,jui+2,:,lb) =  35.*work(iui,jui,:,lb)     - 42.*work(iui,jui-1,:,lb)   + 21.*work(iui,jui-2,:,lb)   -  4.* work(iui,jui-3,:,lb) &
                                - 42.*work(iui-1,jui,:,lb)   + 28.*work(iui-1,jui-1,:,lb) -  6.*work(iui-1,jui-2,:,lb) + 21.* work(iui-2,jui,:,lb) &
                                -  6.*work(iui-2,jui-1,:,lb) -  4.*work(iui-3,jui,:,lb)
      endif
    endif

    if (ndim == 3) then
      if (nbr_blks(5) <= -20) then
        if (extrap /= 0) then
          work(:,:,kli-1,lb) = 2.*work(:,:,kli,lb) - work(:,:,kli+1,lb)
          work(:,:,kli-2,lb) = 3.*work(:,:,kli,lb) - 2.*work(:,:,kli+1,lb)
          work(:,:,kli-3,lb) = 4.*work(:,:,kli,lb) - 3.*work(:,:,kli+1,lb)
          work(:,:,kli-4,lb) = 5.*work(:,:,kli,lb) - 4.*work(:,:,kli+1,lb)
        else
          work(:,:,kli-1,lb) = -work(:,:,kli,lb)
          work(:,:,kli-2,lb) = -3.*work(:,:,kli,lb)
          work(:,:,kli-3,lb) = -5.*work(:,:,kli,lb)
          work(:,:,kli-4,lb) = -7.*work(:,:,kli,lb)
        endif
      endif
      if ((nbr_blks(1) <= -20) .and. (nbr_blks(5) <= -20)) then
        work(ili-1,:,kli-1,lb) = +10*work(ili+0,:,kli+0,lb)-10*work(ili+0,:,kli+1,lb) &
                                  +5*work(ili+0,:,kli+2,lb) -1*work(ili+0,:,kli+3,lb) &
                                 -10*work(ili+1,:,kli+0,lb) +5*work(ili+1,:,kli+1,lb) &
                                  -1*work(ili+1,:,kli+2,lb) +5*work(ili+2,:,kli+0,lb) &
                                  -1*work(ili+2,:,kli+1,lb) -1*work(ili+3,:,kli+0,lb)

        work(ili-1,:,kli-2,lb) = +20*work(ili+0,:,kli+0,lb)-30*work(ili+0,:,kli+1,lb) &
                                 +18*work(ili+0,:,kli+2,lb) -4*work(ili+0,:,kli+3,lb) &
                                 -15*work(ili+1,:,kli+0,lb)+12*work(ili+1,:,kli+1,lb) &
                                  -3*work(ili+1,:,kli+2,lb) +6*work(ili+2,:,kli+0,lb) &
                                  -2*work(ili+2,:,kli+1,lb) -1*work(ili+3,:,kli+0,lb)

        work(ili-2,:,kli-1,lb) = +20*work(ili+0,:,kli+0,lb)-15*work(ili+0,:,kli+1,lb) &
                                  +6*work(ili+0,:,kli+2,lb) -1*work(ili+0,:,kli+3,lb) &
                                 -30*work(ili+1,:,kli+0,lb)+12*work(ili+1,:,kli+1,lb) &
                                  -2*work(ili+1,:,kli+2,lb)+18*work(ili+2,:,kli+0,lb) &
                                  -3*work(ili+2,:,kli+1,lb) -4*work(ili+3,:,kli+0,lb)

        work(ili-2,:,kli-2,lb) = +35*work(ili+0,:,kli+0,lb)-42*work(ili+0,:,kli+1,lb) &
                                 +21*work(ili+0,:,kli+2,lb) -4*work(ili+0,:,kli+3,lb) &
                                 -42*work(ili+1,:,kli+0,lb)+28*work(ili+1,:,kli+1,lb) &
                                  -6*work(ili+1,:,kli+2,lb)+21*work(ili+2,:,kli+0,lb) &
                                  -6*work(ili+2,:,kli+1,lb) -4*work(ili+3,:,kli+0,lb)
      endif
      if ((nbr_blks(2) <= -20) .and. (nbr_blks(5) <= -20)) then
        work(iui+1,:,kli-1,lb) = +10*work(iui-0,:,kli+0,lb)-10*work(iui-0,:,kli+1,lb) &
                                  +5*work(iui-0,:,kli+2,lb) -1*work(iui-0,:,kli+3,lb) &
                                 -10*work(iui-1,:,kli+0,lb) +5*work(iui-1,:,kli+1,lb) &
                                  -1*work(iui-1,:,kli+2,lb) +5*work(iui-2,:,kli+0,lb) &
                                  -1*work(iui-2,:,kli+1,lb) -1*work(iui-3,:,kli+0,lb)

        work(iui+1,:,kli-2,lb) = +20*work(iui-0,:,kli+0,lb)-30*work(iui-0,:,kli+1,lb) &
                                 +18*work(iui-0,:,kli+2,lb) -4*work(iui-0,:,kli+3,lb) &
                                 -15*work(iui-1,:,kli+0,lb)+12*work(iui-1,:,kli+1,lb) &
                                  -3*work(iui-1,:,kli+2,lb) +6*work(iui-2,:,kli+0,lb) &
                                  -2*work(iui-2,:,kli+1,lb) -1*work(iui-3,:,kli+0,lb)

        work(iui+2,:,kli-1,lb) = +20*work(iui-0,:,kli+0,lb)-15*work(iui-0,:,kli+1,lb) &
                                  +6*work(iui-0,:,kli+2,lb) -1*work(iui-0,:,kli+3,lb) &
                                 -30*work(iui-1,:,kli+0,lb)+12*work(iui-1,:,kli+1,lb) &
                                  -2*work(iui-1,:,kli+2,lb)+18*work(iui-2,:,kli+0,lb) &
                                  -3*work(iui-2,:,kli+1,lb) -4*work(iui-3,:,kli+0,lb)

        work(iui+2,:,kli-2,lb) = +35*work(iui-0,:,kli+0,lb)-42*work(iui-0,:,kli+1,lb) &
                                 +21*work(iui-0,:,kli+2,lb) -4*work(iui-0,:,kli+3,lb) &
                                 -42*work(iui-1,:,kli+0,lb)+28*work(iui-1,:,kli+1,lb) &
                                  -6*work(iui-1,:,kli+2,lb)+21*work(iui-2,:,kli+0,lb) &
                                  -6*work(iui-2,:,kli+1,lb) -4*work(iui-3,:,kli+0,lb)
      endif
      if ((nbr_blks(3) <= -20) .and. (nbr_blks(5) <= -20)) then
        work(:,jli-1,kli-1,lb) = +10*work(:,jli+0,kli+0,lb)-10*work(:,jli+0,kli+1,lb) &
                                  +5*work(:,jli+0,kli+2,lb) -1*work(:,jli+0,kli+3,lb) &
                                 -10*work(:,jli+1,kli+0,lb) +5*work(:,jli+1,kli+1,lb) &
                                  -1*work(:,jli+1,kli+2,lb) +5*work(:,jli+2,kli+0,lb) &
                                  -1*work(:,jli+2,kli+1,lb) -1*work(:,jli+3,kli+0,lb)

        work(:,jli-1,kli-2,lb) = +20*work(:,jli+0,kli+0,lb)-30*work(:,jli+0,kli+1,lb) &
                                 +18*work(:,jli+0,kli+2,lb) -4*work(:,jli+0,kli+3,lb) &
                                 -15*work(:,jli+1,kli+0,lb)+12*work(:,jli+1,kli+1,lb) &
                                  -3*work(:,jli+1,kli+2,lb) +6*work(:,jli+2,kli+0,lb) &
                                  -2*work(:,jli+2,kli+1,lb) -1*work(:,jli+3,kli+0,lb)

        work(:,jli-2,kli-1,lb) = +20*work(:,jli+0,kli+0,lb)-15*work(:,jli+0,kli+1,lb) &
                                  +6*work(:,jli+0,kli+2,lb) -1*work(:,jli+0,kli+3,lb) &
                                 -30*work(:,jli+1,kli+0,lb)+12*work(:,jli+1,kli+1,lb) &
                                  -2*work(:,jli+1,kli+2,lb)+18*work(:,jli+2,kli+0,lb) &
                                  -3*work(:,jli+2,kli+1,lb) -4*work(:,jli+3,kli+0,lb)

        work(:,jli-2,kli-2,lb) = +35*work(:,jli+0,kli+0,lb)-42*work(:,jli+0,kli+1,lb) &
                                 +21*work(:,jli+0,kli+2,lb) -4*work(:,jli+0,kli+3,lb) &
                                 -42*work(:,jli+1,kli+0,lb)+28*work(:,jli+1,kli+1,lb) &
                                  -6*work(:,jli+1,kli+2,lb)+21*work(:,jli+2,kli+0,lb) &
                                  -6*work(:,jli+2,kli+1,lb) -4*work(:,jli+3,kli+0,lb)
      endif
      if ((nbr_blks(4) <= -20) .and. (nbr_blks(5) <= -20)) then
        work(:,jui+1,kli-1,lb) = +10*work(:,jui-0,kli+0,lb)-10*work(:,jui-0,kli+1,lb) &
                                  +5*work(:,jui-0,kli+2,lb) -1*work(:,jui-0,kli+3,lb) &
                                 -10*work(:,jui-1,kli+0,lb) +5*work(:,jui-1,kli+1,lb) &
                                  -1*work(:,jui-1,kli+2,lb) +5*work(:,jui-2,kli+0,lb) &
                                  -1*work(:,jui-2,kli+1,lb) -1*work(:,jui-3,kli+0,lb)

        work(:,jui+1,kli-2,lb) = +20*work(:,jui-0,kli+0,lb)-30*work(:,jui-0,kli+1,lb) &
                                 +18*work(:,jui-0,kli+2,lb) -4*work(:,jui-0,kli+3,lb) &
                                 -15*work(:,jui-1,kli+0,lb)+12*work(:,jui-1,kli+1,lb) &
                                  -3*work(:,jui-1,kli+2,lb) +6*work(:,jui-2,kli+0,lb) &
                                  -2*work(:,jui-2,kli+1,lb) -1*work(:,jui-3,kli+0,lb)

        work(:,jui+2,kli-1,lb) = +20*work(:,jui-0,kli+0,lb)-15*work(:,jui-0,kli+1,lb) &
                                  +6*work(:,jui-0,kli+2,lb) -1*work(:,jui-0,kli+3,lb) &
                                 -30*work(:,jui-1,kli+0,lb)+12*work(:,jui-1,kli+1,lb) &
                                  -2*work(:,jui-1,kli+2,lb)+18*work(:,jui-2,kli+0,lb) &
                                  -3*work(:,jui-2,kli+1,lb) -4*work(:,jui-3,kli+0,lb)

        work(:,jui+2,kli-2,lb) = +35*work(:,jui-0,kli+0,lb)-42*work(:,jui-0,kli+1,lb) &
                                 +21*work(:,jui-0,kli+2,lb) -4*work(:,jui-0,kli+3,lb) &
                                 -42*work(:,jui-1,kli+0,lb)+28*work(:,jui-1,kli+1,lb) &
                                  -6*work(:,jui-1,kli+2,lb)+21*work(:,jui-2,kli+0,lb) &
                                  -6*work(:,jui-2,kli+1,lb) -4*work(:,jui-3,kli+0,lb)
      endif

! Should use the correct 2nd-order expressions for external corners, but for
! now we are just setting these cells to zero.

      if ((nbr_blks(1) <= -20) .and. (nbr_blks(3) <= -20) .and. (nbr_blks(5) <= -20)) then
        work(ili-1,jli-1,kli-1,lb) = 0.
        work(ili-2,jli-1,kli-1,lb) = 0.
        work(ili-1,jli-2,kli-1,lb) = 0.
        work(ili-2,jli-2,kli-1,lb) = 0.
        work(ili-1,jli-1,kli-2,lb) = 0.
        work(ili-2,jli-1,kli-2,lb) = 0.
        work(ili-1,jli-2,kli-2,lb) = 0.
        work(ili-2,jli-2,kli-2,lb) = 0.
      endif
      if ((nbr_blks(1) <= -20) .and. (nbr_blks(4) <= -20) .and. (nbr_blks(5) <= -20)) then
        work(ili-1,jui+1,kli-1,lb) = 0.
        work(ili-2,jui+1,kli-1,lb) = 0.
        work(ili-1,jui+2,kli-1,lb) = 0.
        work(ili-2,jui+2,kli-1,lb) = 0.
        work(ili-1,jui+1,kli-2,lb) = 0.
        work(ili-2,jui+1,kli-2,lb) = 0.
        work(ili-1,jui+2,kli-2,lb) = 0.
        work(ili-2,jui+2,kli-2,lb) = 0.
      endif
      if ((nbr_blks(2) <= -20) .and. (nbr_blks(3) <= -20) .and. (nbr_blks(5) <= -20)) then
        work(iui+1,jli-1,kli-1,lb) = 0.
        work(iui+2,jli-1,kli-1,lb) = 0.
        work(iui+1,jli-2,kli-1,lb) = 0.
        work(iui+2,jli-2,kli-1,lb) = 0.
        work(iui+1,jli-1,kli-2,lb) = 0.
        work(iui+2,jli-1,kli-2,lb) = 0.
        work(iui+1,jli-2,kli-2,lb) = 0.
        work(iui+2,jli-2,kli-2,lb) = 0.
      endif
      if ((nbr_blks(2) <= -20) .and. (nbr_blks(4) <= -20) .and. (nbr_blks(5) <= -20)) then
        work(iui+1,jui+1,kli-1,lb) = 0.
        work(iui+2,jui+1,kli-1,lb) = 0.
        work(iui+1,jui+2,kli-1,lb) = 0.
        work(iui+2,jui+2,kli-1,lb) = 0.
        work(iui+1,jui+1,kli-2,lb) = 0.
        work(iui+2,jui+1,kli-2,lb) = 0.
        work(iui+1,jui+2,kli-2,lb) = 0.
        work(iui+2,jui+2,kli-2,lb) = 0.
      endif


      if (nbr_blks(6) <= -20) then
        if (extrap /= 0) then
          work(:,:,kui+1,lb) = 2.*work(:,:,kui,lb) - work(:,:,kui-1,lb)
          work(:,:,kui+2,lb) = 3.*work(:,:,kui,lb) - 2.*work(:,:,kui-1,lb)
          work(:,:,kui+3,lb) = 4.*work(:,:,kui,lb) - 3.*work(:,:,kui-1,lb)
          work(:,:,kui+4,lb) = 5.*work(:,:,kui,lb) - 4.*work(:,:,kui-1,lb)
        else
          work(:,:,kui+1,lb) = -work(:,:,kui,lb)
          work(:,:,kui+2,lb) = -3.*work(:,:,kui,lb)
          work(:,:,kui+3,lb) = -5.*work(:,:,kui,lb)
          work(:,:,kui+4,lb) = -7.*work(:,:,kui,lb)
        endif
      endif
      if ((nbr_blks(1) <= -20) .and. (nbr_blks(6) <= -20)) then
        work(ili-1,:,kui+1,lb) = +10*work(ili+0,:,kui-0,lb)-10*work(ili+0,:,kui-1,lb) &
                                  +5*work(ili+0,:,kui-2,lb) -1*work(ili+0,:,kui-3,lb) &
                                 -10*work(ili+1,:,kui-0,lb) +5*work(ili+1,:,kui-1,lb) &
                                  -1*work(ili+1,:,kui-2,lb) +5*work(ili+2,:,kui-0,lb) &
                                  -1*work(ili+2,:,kui-1,lb) -1*work(ili+3,:,kui-0,lb)

        work(ili-1,:,kui+2,lb) = +20*work(ili+0,:,kui-0,lb)-30*work(ili+0,:,kui-1,lb) &
                                 +18*work(ili+0,:,kui-2,lb) -4*work(ili+0,:,kui-3,lb) &
                                 -15*work(ili+1,:,kui-0,lb)+12*work(ili+1,:,kui-1,lb) &
                                  -3*work(ili+1,:,kui-2,lb) +6*work(ili+2,:,kui-0,lb) &
                                  -2*work(ili+2,:,kui-1,lb) -1*work(ili+3,:,kui-0,lb)

        work(ili-2,:,kui+1,lb) = +20*work(ili+0,:,kui-0,lb)-15*work(ili+0,:,kui-1,lb) &
                                  +6*work(ili+0,:,kui-2,lb) -1*work(ili+0,:,kui-3,lb) &
                                 -30*work(ili+1,:,kui-0,lb)+12*work(ili+1,:,kui-1,lb) &
                                  -2*work(ili+1,:,kui-2,lb)+18*work(ili+2,:,kui-0,lb) &
                                  -3*work(ili+2,:,kui-1,lb) -4*work(ili+3,:,kui-0,lb)

        work(ili-2,:,kui+2,lb) = +35*work(ili+0,:,kui-0,lb)-42*work(ili+0,:,kui-1,lb) &
                                 +21*work(ili+0,:,kui-2,lb) -4*work(ili+0,:,kui-3,lb) &
                                 -42*work(ili+1,:,kui-0,lb)+28*work(ili+1,:,kui-1,lb) &
                                  -6*work(ili+1,:,kui-2,lb)+21*work(ili+2,:,kui-0,lb) &
                                  -6*work(ili+2,:,kui-1,lb) -4*work(ili+3,:,kui-0,lb)
      endif
      if ((nbr_blks(2) <= -20) .and. (nbr_blks(6) <= -20)) then
        work(iui+1,:,kui+1,lb) = +10*work(iui-0,:,kui-0,lb)-10*work(iui-0,:,kui-1,lb) &
                                  +5*work(iui-0,:,kui-2,lb) -1*work(iui-0,:,kui-3,lb) &
                                 -10*work(iui-1,:,kui-0,lb) +5*work(iui-1,:,kui-1,lb) &
                                  -1*work(iui-1,:,kui-2,lb) +5*work(iui-2,:,kui-0,lb) &
                                  -1*work(iui-2,:,kui-1,lb) -1*work(iui-3,:,kui-0,lb)

        work(iui+1,:,kui+2,lb) = +20*work(iui-0,:,kui-0,lb)-30*work(iui-0,:,kui-1,lb) &
                                 +18*work(iui-0,:,kui-2,lb) -4*work(iui-0,:,kui-3,lb) &
                                 -15*work(iui-1,:,kui-0,lb)+12*work(iui-1,:,kui-1,lb) &
                                  -3*work(iui-1,:,kui-2,lb) +6*work(iui-2,:,kui-0,lb) &
                                  -2*work(iui-2,:,kui-1,lb) -1*work(iui-3,:,kui-0,lb)

        work(iui+2,:,kui+1,lb) = +20*work(iui-0,:,kui-0,lb)-15*work(iui-0,:,kui-1,lb) &
                                  +6*work(iui-0,:,kui-2,lb) -1*work(iui-0,:,kui-3,lb) &
                                 -30*work(iui-1,:,kui-0,lb)+12*work(iui-1,:,kui-1,lb) &
                                  -2*work(iui-1,:,kui-2,lb)+18*work(iui-2,:,kui-0,lb) &
                                  -3*work(iui-2,:,kui-1,lb) -4*work(iui-3,:,kui-0,lb)

        work(iui+2,:,kui+2,lb) = +35*work(iui-0,:,kui-0,lb)-42*work(iui-0,:,kui-1,lb) &
                                 +21*work(iui-0,:,kui-2,lb) -4*work(iui-0,:,kui-3,lb) &
                                 -42*work(iui-1,:,kui-0,lb)+28*work(iui-1,:,kui-1,lb) &
                                  -6*work(iui-1,:,kui-2,lb)+21*work(iui-2,:,kui-0,lb) &
                                  -6*work(iui-2,:,kui-1,lb) -4*work(iui-3,:,kui-0,lb)
      endif
      if ((nbr_blks(3) <= -20) .and. (nbr_blks(6) <= -20)) then
        work(:,jli-1,kui+1,lb) = +10*work(:,jli+0,kui-0,lb)-10*work(:,jli+0,kui-1,lb) &
                                  +5*work(:,jli+0,kui-2,lb) -1*work(:,jli+0,kui-3,lb) &
                                 -10*work(:,jli+1,kui-0,lb) +5*work(:,jli+1,kui-1,lb) &
                                  -1*work(:,jli+1,kui-2,lb) +5*work(:,jli+2,kui-0,lb) &
                                  -1*work(:,jli+2,kui-1,lb) -1*work(:,jli+3,kui-0,lb)

        work(:,jli-1,kui+2,lb) = +20*work(:,jli+0,kui-0,lb)-30*work(:,jli+0,kui-1,lb) &
                                 +18*work(:,jli+0,kui-2,lb) -4*work(:,jli+0,kui-3,lb) &
                                 -15*work(:,jli+1,kui-0,lb)+12*work(:,jli+1,kui-1,lb) &
                                  -3*work(:,jli+1,kui-2,lb) +6*work(:,jli+2,kui-0,lb) &
                                  -2*work(:,jli+2,kui-1,lb) -1*work(:,jli+3,kui-0,lb)

        work(:,jli-2,kui+1,lb) = +20*work(:,jli+0,kui-0,lb)-15*work(:,jli+0,kui-1,lb) &
                                  +6*work(:,jli+0,kui-2,lb) -1*work(:,jli+0,kui-3,lb) &
                                 -30*work(:,jli+1,kui-0,lb)+12*work(:,jli+1,kui-1,lb) &
                                  -2*work(:,jli+1,kui-2,lb)+18*work(:,jli+2,kui-0,lb) &
                                  -3*work(:,jli+2,kui-1,lb) -4*work(:,jli+3,kui-0,lb)

        work(:,jli-2,kui+2,lb) = +35*work(:,jli+0,kui-0,lb)-42*work(:,jli+0,kui-1,lb) &
                                 +21*work(:,jli+0,kui-2,lb) -4*work(:,jli+0,kui-3,lb) &
                                 -42*work(:,jli+1,kui-0,lb)+28*work(:,jli+1,kui-1,lb) &
                                  -6*work(:,jli+1,kui-2,lb)+21*work(:,jli+2,kui-0,lb) &
                                  -6*work(:,jli+2,kui-1,lb) -4*work(:,jli+3,kui-0,lb)
      endif
      if ((nbr_blks(4) <= -20) .and. (nbr_blks(6) <= -20)) then
        work(:,jui+1,kui+1,lb) = +10*work(:,jui-0,kui-0,lb)-10*work(:,jui-0,kui-1,lb) &
                                  +5*work(:,jui-0,kui-2,lb) -1*work(:,jui-0,kui-3,lb) &
                                 -10*work(:,jui-1,kui-0,lb) +5*work(:,jui-1,kui-1,lb) &
                                  -1*work(:,jui-1,kui-2,lb) +5*work(:,jui-2,kui-0,lb) &
                                  -1*work(:,jui-2,kui-1,lb) -1*work(:,jui-3,kui-0,lb)

        work(:,jui+1,kui+2,lb) = +20*work(:,jui-0,kui-0,lb)-30*work(:,jui-0,kui-1,lb) &
                                 +18*work(:,jui-0,kui-2,lb) -4*work(:,jui-0,kui-3,lb) &
                                 -15*work(:,jui-1,kui-0,lb)+12*work(:,jui-1,kui-1,lb) &
                                  -3*work(:,jui-1,kui-2,lb) +6*work(:,jui-2,kui-0,lb) &
                                  -2*work(:,jui-2,kui-1,lb) -1*work(:,jui-3,kui-0,lb)

        work(:,jui+2,kui+1,lb) = +20*work(:,jui-0,kui-0,lb)-15*work(:,jui-0,kui-1,lb) &
                                  +6*work(:,jui-0,kui-2,lb) -1*work(:,jui-0,kui-3,lb) &
                                 -30*work(:,jui-1,kui-0,lb)+12*work(:,jui-1,kui-1,lb) &
                                  -2*work(:,jui-1,kui-2,lb)+18*work(:,jui-2,kui-0,lb) &
                                  -3*work(:,jui-2,kui-1,lb) -4*work(:,jui-3,kui-0,lb)

        work(:,jui+2,kui+2,lb) = +35*work(:,jui-0,kui-0,lb)-42*work(:,jui-0,kui-1,lb) &
                                 +21*work(:,jui-0,kui-2,lb) -4*work(:,jui-0,kui-3,lb) &
                                 -42*work(:,jui-1,kui-0,lb)+28*work(:,jui-1,kui-1,lb) &
                                  -6*work(:,jui-1,kui-2,lb)+21*work(:,jui-2,kui-0,lb) &
                                  -6*work(:,jui-2,kui-1,lb) -4*work(:,jui-3,kui-0,lb)
      endif

! Should use the correct 2nd-order expressions for external corners, but for
! now we are just setting these cells to zero.

      if ((nbr_blks(1) <= -20) .and. (nbr_blks(3) <= -20) .and. (nbr_blks(6) <= -20)) then
        work(ili-1,jli-1,kui+1,lb) = 0.
        work(ili-2,jli-1,kui+1,lb) = 0.
        work(ili-1,jli-2,kui+1,lb) = 0.
        work(ili-2,jli-2,kui+1,lb) = 0.
        work(ili-1,jli-1,kui+2,lb) = 0.
        work(ili-2,jli-1,kui+2,lb) = 0.
        work(ili-1,jli-2,kui+2,lb) = 0.
        work(ili-2,jli-2,kui+2,lb) = 0.
      endif
      if ((nbr_blks(1) <= -20) .and. (nbr_blks(4) <= -20) .and. (nbr_blks(6) <= -20)) then
        work(ili-1,jui+1,kui+1,lb) = 0.
        work(ili-2,jui+1,kui+1,lb) = 0.
        work(ili-1,jui+2,kui+1,lb) = 0.
        work(ili-2,jui+2,kui+1,lb) = 0.
        work(ili-1,jui+1,kui+2,lb) = 0.
        work(ili-2,jui+1,kui+2,lb) = 0.
        work(ili-1,jui+2,kui+2,lb) = 0.
        work(ili-2,jui+2,kui+2,lb) = 0.
      endif
      if ((nbr_blks(2) <= -20) .and. (nbr_blks(3) <= -20) .and. (nbr_blks(6) <= -20)) then
        work(iui+1,jli-1,kui+1,lb) = 0.
        work(iui+2,jli-1,kui+1,lb) = 0.
        work(iui+1,jli-2,kui+1,lb) = 0.
        work(iui+2,jli-2,kui+1,lb) = 0.
        work(iui+1,jli-1,kui+2,lb) = 0.
        work(iui+2,jli-1,kui+2,lb) = 0.
        work(iui+1,jli-2,kui+2,lb) = 0.
        work(iui+2,jli-2,kui+2,lb) = 0.
      endif
      if ((nbr_blks(2) <= -20) .and. (nbr_blks(4) <= -20) .and. (nbr_blks(6) <= -20)) then
        work(iui+1,jui+1,kui+1,lb) = 0.
        work(iui+2,jui+1,kui+1,lb) = 0.
        work(iui+1,jui+2,kui+1,lb) = 0.
        work(iui+2,jui+2,kui+1,lb) = 0.
        work(iui+1,jui+1,kui+2,lb) = 0.
        work(iui+2,jui+1,kui+2,lb) = 0.
        work(iui+1,jui+2,kui+2,lb) = 0.
        work(iui+2,jui+2,kui+2,lb) = 0.
      endif
    endif

  enddo

else if (mg_bnd_cond == MG_BND_NEUMANN) then ! Neumann

  do lb = 1, lnblocks

    nbr_blks = dBaseNeighborBlockList(lb)

    if (nbr_blks(1) <= -20) work(ili-1,:,:,lb) = work(ili,:,:,lb)
    if (nbr_blks(2) <= -20) work(iui+1,:,:,lb) = work(iui,:,:,lb)

    if (ndim >= 2) then
      if (nbr_blks(3) <= -20) work(:,jli-1,:,lb) = work(:,jli,:,lb)
      if (nbr_blks(4) <= -20) work(:,jui+1,:,lb) = work(:,jui,:,lb)
    endif

    if (ndim == 3) then
      if (nbr_blks(5) <= -20) work(:,:,kli-1,lb) = work(:,:,kli,lb)
      if (nbr_blks(6) <= -20) work(:,:,kui+1,lb) = work(:,:,kui,lb)
    endif

  enddo

endif

!===============================================================================

return
end 

