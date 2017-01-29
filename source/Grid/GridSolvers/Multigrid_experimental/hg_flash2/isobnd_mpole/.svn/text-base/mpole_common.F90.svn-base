!!****if* source/Grid/GridSolvers/Multigrid_experimental/hg_flash2/isobnd_mpole/mpole_common
!!
!! NAME
!!
!!  mpole_common
!!
!!
!! SYNOPSIS
!!
!!  use mpole_common
!!
!!
!! DESCRIPTION
!!
!!  Variable declarations for the multipole Poisson solver.
!!
!!  Multipole module data structures -- these are allocated in init_mpole 
!!  once the size of the mesh is computed.  These are then used to share 
!!  data across the various mpole functions
!!
!!    Moment(q,k,p,l,m)     The (l,m)-th moment (k=1, even; 2, odd;
!!                          p=1, inner; 2, outer) of the image
!!                          mass distribution as a fctn of q*dsinv
!!
!!    Momtmp(q)             Temporary array to receive radial
!!                          samples of Moment() in parallel
!!                          reduction operation
!!
!!    dsinv                 Moment array sample spacing
!!
!!    mpole_lmax            Maximum multipole moment (runtime set)
!!
!!    mpole_mmax            Maximum azimuthal moment (determined by geometry)
!!
!!    qmax                  Maximum number of radial samples
!!
!!    Mtot                  Total mass
!!
!!    X/Y/Zcm               Location of center of mass
!!
!!    mpole_geometry        Geometry type for grid
!!
!!    costable(m)           Table containing the cosine of
!!                          m * the current azimuthal angle; also
!!                          sintable(m)
!!
!!    rpower(l)             The current cell's density*rprime^l
!!
!!    rprinv(l)             The current cell's density*rprime^-(l+1)
!!
!!    Leg_fact(l,m)         Factorial normalization coefficients
!!                          for the associated Legendre function
!!
!!***

module mpole_common

  implicit none
  save


! database keys -- these will be filled by init_mpole from the FLASH
! database.

  integer :: ndim

  integer :: izn, iznl, iznr
  integer :: ixCoord, iyCoord, izCoord

  integer :: MyPE, MasterPe

  real    :: xmin, xmax, ymin, ymax, zmin, zmax

  integer :: lrefine_max, Nblockx, Nblocky, Nblockz

  integer :: nxb, nyb, nzb, imin, imax, jmin, jmax, kmin, kmax, k2d, k3d


! Physical/mathematical constants -- these are filling in init_mpole

  real    :: pi, twopi, fourpi


  real               :: Mtot, Xcm, Ycm, Zcm, dsinv

  real, allocatable  :: Moment(:,:,:,:,:), Momtmp(:)

  real, allocatable  :: costable(:), sintable(:)

  real, allocatable  :: rpower(:), rprinv(:)

  real, allocatable  :: Legk1(:,:), Legk2(:,:), Leg_fact(:,:)

  integer            :: qmax, mpole_lmax, mpole_mmax, mpole_geometry


! Constants used to index moment array

  integer, parameter :: Even = 1, Odd = 2, Inner = 1, Outer = 2


! Supported geometry constants -- we do a geometry check to make sure
! that we support whatever the request geometry is in init_mpole

  integer, parameter :: G_3DCARTESIAN  = 1, G_2DCYLINDRICAL = 2, &
                        G_1DSPHERICAL  = 3, G_3DSPHERICAL   = 4, &
                        G_3DCYLINDRICAL= 5


! in 2-d cylindrical coordinates, we allow a single quadrant of the star
! to be simulated if quadrant = .true.

  logical :: quadrant


end module mpole_common
