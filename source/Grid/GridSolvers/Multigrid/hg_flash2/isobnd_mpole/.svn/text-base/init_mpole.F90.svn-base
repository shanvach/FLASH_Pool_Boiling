!!****if* source/Grid/GridSolvers/Multigrid/hg_flash2/isobnd_mpole/init_mpole
!!
!! NAME
!!
!!  init_mpole
!!
!! 
!! SYNOPSIS
!!
!!  init_mpole()
!!
!!
!! DESCRIPTION
!!
!!  Initialize the multipole Poisson solver.  Read in any of the
!!  runtime parameters for this solver.  All solver common data
!!  is stored in the mpole_common module
!!
!!***

subroutine init_mpole

  use dBase, ONLY: dBaseKeyNumber, dBasePropertyInteger, dBasePropertyReal, &
                   CARTESIAN, CYLINDRICAL, SPHERICAL

  use mpole_common

  use runtime_parameters

  use physical_constants

  use logfile, ONLY: write_logfile

  implicit none

  real    :: Nxmax, Nymax, Nzmax, dxmin, dymin, dzmin
  real    :: factrl
  integer :: istat, m, l

  integer :: meshGeom

  character(len=128) :: str_buffer
!===============================================================================

!               Initialize database keys, constants, etc.

  MyPE     = dBasePropertyInteger("MyProcessor")
  MasterPE = dBasePropertyInteger("MasterProcessor")

  ndim = dBasePropertyInteger("Dimensionality")

  k2d = 0
  if (ndim >= 2) k2d = 1
  k3d = 0
  if (ndim == 3) k3d = 1
  
  nxb  = dBasePropertyInteger("xBlockSize")
  nyb  = dBasePropertyInteger("yBlockSize")
  nzb  = dBasePropertyInteger("zBlockSize")

  imin = dBasePropertyInteger("xBlockLo")
  imax = dBasePropertyInteger("xBlockHi")
  jmin = dBasePropertyInteger("yBlockLo")
  jmax = dBasePropertyInteger("yBlockHi")
  kmin = dBasePropertyInteger("zBlockLo")
  kmax = dBasePropertyInteger("zBlockHi")

  call get_constant_from_db ("pi", pi)
  twopi  = 2.*pi
  fourpi = 4.*pi

  call get_parm_from_context(global_parm_context, "mpole_lmax", mpole_lmax)

  call get_parm_from_context(global_parm_context, "quadrant", quadrant)

  call get_parm_from_context(global_parm_context, "lrefine_max", lrefine_max)

  call get_parm_from_context(global_parm_context, "Nblockx", Nblockx)
  call get_parm_from_context(global_parm_context, "Nblocky", Nblocky)
  call get_parm_from_context(global_parm_context, "Nblockz", Nblockz)

  call get_parm_from_context(global_parm_context, "xmin", xmin)
  call get_parm_from_context(global_parm_context, "xmax", xmax)
  call get_parm_from_context(global_parm_context, "ymin", ymin)
  call get_parm_from_context(global_parm_context, "ymax", ymax)
  call get_parm_from_context(global_parm_context, "zmin", zmin)
  call get_parm_from_context(global_parm_context, "zmax", zmax)

  izn   = dBaseKeyNumber("zn")
  iznl  = dBaseKeyNumber("znl")
  iznr  = dBaseKeyNumber("znr")

  ixCoord = dBaseKeyNumber("xCoord")
  iyCoord = dBaseKeyNumber("yCoord")
  izCoord = dBaseKeyNumber("zCoord")

  meshGeom = dBasePropertyInteger("MeshGeometry")

!               Check if we support the requested grid geometry.
!               Only bounded grid geometries are supported, and
!               support for the two which use angular coordinates
!               (3D cylindrical and 3D spherical) is deferred for now.

  if ((ndim == 3) .and. (meshGeom == CARTESIAN)) then

     mpole_geometry = G_3DCARTESIAN
     mpole_mmax = mpole_lmax
 
  else if ((ndim == 2) .and. (meshGeom == CYLINDRICAL)) then

     mpole_geometry = G_2DCYLINDRICAL
     mpole_mmax = 0
     call write_logfile ('init_mpole:  2d axisymmetry, ignoring m > 0 moments')

  else if ((ndim == 1) .and. (meshGeom == SPHERICAL)) then

     mpole_geometry = G_1DSPHERICAL
     mpole_lmax = 0
     mpole_mmax = 0
     call write_logfile &
       ('init_mpole:  1d spherical symmetry, ignoring l > 0 moments')

  else

     call abort_flash ('init_mpole:  FATAL:  unsupported geometry!')

  endif

! we are only allowed to do a quadrant (i.e. enforce reflection symmetry
! about y=0) if we are in 2-d cylindrical coords
  if ((quadrant) .AND. (mpole_geometry /= G_2DCYLINDRICAL)) then
     call abort_flash('ERROR: quadrant only allowed in 2-d cylindrical geometry')
  endif
    

!===============================================================================

!               Maximum number of zones across each dimension, if each
!               dimension were to become fully refined.  Also compute
!               minimum zone spacings at the maximum refinement level.
!               Assume the domain is a Cartesian box.

  Nxmax = Nblockx * nxb * 2.**(lrefine_max-1)

  dxmin = (xmax - xmin) / Nxmax

  if (ndim >= 2) then
     Nymax = Nblocky * nyb * 2.**(lrefine_max-1)
     dymin = (ymax - ymin) / Nymax
  else
     Nymax = 0
     dymin = 1.
  endif

  if (ndim == 3) then
     Nzmax = Nblockz * nzb * 2.**(lrefine_max-1)
     dzmin = (zmax - zmin) / Nzmax
  else
     Nzmax = 0
     dzmin = 1.
  endif

!               Inverse of sample spacing to use for moment arrays.

  dsinv = 2. / (dxmin*dymin*dzmin)**(1./ndim)

!               Number of radial samples to use in moment arrays.

  qmax = 2 * int( sqrt(Nxmax**2 + Nymax**2 + Nzmax**2) ) + 3

  if (MyPE == MasterPE) then
     write (str_buffer,*) 'init_mpole:  using ', qmax, ' radial samples', &
       & ' moment array:  ', & 
       & qmax*2*2*(mpole_lmax+1)*(mpole_mmax+1), ' items'
     call write_logfile(str_buffer)
  endif

!               Allocate moment arrays and other data structures.

  allocate ( Moment(0:qmax,1:2,1:2,0:mpole_lmax,0:mpole_mmax), stat=istat)
  if (istat>0) call abort_flash ("init_mpole: Moment() allocate failed!")

  allocate( Momtmp(0:qmax), stat=istat)
  if (istat>0) call abort_flash ("init_mpole: Momtmp() allocate failed!")

  allocate( costable(0:mpole_mmax), stat=istat)
  if (istat>0) call abort_flash ("init_mpole: costable() allocate failed!")

  allocate( sintable(0:mpole_mmax), stat=istat)
  if (istat>0) call abort_flash ("init_mpole: sintable() allocate failed!")

  allocate( rpower(0:mpole_lmax), stat=istat)
  if (istat>0) call abort_flash ("init_mpole: rpower() allocate failed!")

  allocate( Legk1(0:mpole_lmax,0:mpole_mmax), stat=istat)
  if (istat>0) call abort_flash ("init_mpole: Legk1() allocate failed!")

  allocate( rprinv(0:mpole_lmax), stat=istat)
  if (istat>0) call abort_flash ("init_mpole: rprinv() allocate failed!")

  allocate( Legk2(0:mpole_lmax,0:mpole_mmax), stat=istat)
  if (istat>0) call abort_flash ("init_mpole: Legk2() allocate failed!")

  allocate( Leg_fact(0:mpole_lmax,0:mpole_mmax), stat=istat)
  if (istat>0) call abort_flash ("init_mpole: leg_fact() allocate failed!")

!                       Coefficients for Legendre polynomials.

  do m = 0, mpole_mmax
    do l = m+2, mpole_lmax
      Legk1(l,m) = real(2*l - 1) / real(l - m)
      Legk2(l,m) = real(l + m - 1) / real(l - m)
    enddo
  enddo

  if ((mpole_lmax > 0) .and. (mpole_mmax > 0)) then

    do l = 1, mpole_lmax
      factrl = 2.
      do m = 1, l
        factrl = factrl / real((l+m) * (l-m+1))
        Leg_fact(l,m) = factrl
      enddo
    enddo

  endif

!===============================================================================

  return
end subroutine init_mpole
