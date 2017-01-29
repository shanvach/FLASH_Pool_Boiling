!!****if* source/Simulation/SimulationMain/ShockCyl/Simulation_init
!!
!! NAME
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!  Simulation_init()
!!
!! DESCRIPTION   
!!     Initialize all the parameters needed for the Shock Cylinder simulation
!!
!! ARGUMENTS
!!      None.  All data passed through Simulation_data
!!
!! PARAMETERS   
!!      Described in the Config file
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp
  use Multispecies_interface, ONLY : Multispecies_getProperty
  use Grid_interface, ONLY : Grid_getNumProcs
  use sim_ranluxModule

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"
#include "Eos.h"

  
  integer :: i, j, status
  ! for randomization that Cal addeded recently
  integer          :: n_pert
  real             :: d_pert, sdim_i
  integer          :: numPEs

  call Driver_getMype(MESH_COMM, sim_meshMe)

  !! These are only used in Paramesh
#ifdef FLASH_GRID_PARAMESH
  call RuntimeParameters_get("nrefs", nrefs) 
#endif
  call RuntimeParameters_get( 'mach',    mach)

  call RuntimeParameters_get( 'xctr',    xctr )
  call RuntimeParameters_get( 'yctr',    yctr )

  call RuntimeParameters_get( 'ref_rect_x', ref_rect_x)
  call RuntimeParameters_get( 'ref_rect_y', ref_rect_y)

  call RuntimeParameters_get( 'xmax',    xmax )
  call RuntimeParameters_get( 'xmin',    xmin )
  call RuntimeParameters_get( 'ymax',    ymax )
  call RuntimeParameters_get( 'ymin',    ymin )
  if(NDIM == 3) then
     call RuntimeParameters_get( 'zmax',    zmax )
     call RuntimeParameters_get( 'zmin',    zmin )
  endif

  call RuntimeParameters_get("p_amb",   p_amb )
  call RuntimeParameters_get("rho_amb", rho_amb )
  call RuntimeParameters_get("vx_amb",  vx_amb )
  call RuntimeParameters_get("mach",    mach )

  call RuntimeParameters_get( 'sim_xShock',  sim_xShock )

  call RuntimeParameters_get( "vz_sf6",  vz_sf6 )
  call RuntimeParameters_get( 'maxconc', maxconc )

  call RuntimeParameters_get( 'sim_useRadialFit', sim_useRadialFit )
  if(sim_useRadialFit) then
     call RuntimeParameters_get( 'sim_radialFitRadius',  sim_radialFitRadius )
  endif


  call RuntimeParameters_get( 'sim_useRawData', sim_useRawData )
  if (sim_useRawData) then
     call RuntimeParameters_get( 'sim_rawNumPixelsX',   sim_rawNumPixelsX )
     call RuntimeParameters_get( 'sim_rawNumPixelsy',   sim_rawNumPixelsY )
     call RuntimeParameters_get( 'sim_rawMinX', sim_rawMinX )
     call RuntimeParameters_get( 'sim_rawMinY', sim_rawMinY )
  endif

  call RuntimeParameters_get( 'use_rz_sim_data', use_rz_sim_data )
  if (use_rz_sim_data) then

     call RuntimeParameters_get( 'rz_rmax', rz_rmax)
     call RuntimeParameters_get( 'rz_zmax', rz_zmax)
     call RuntimeParameters_get( 'rz_3d_use_sym', rz_3d_use_sym )
     call RuntimeParameters_get( 'rz_pert_amp',   rz_pert_amp )
     call RuntimeParameters_get( 'rz_pert_zlen',  rz_pert_zlen )

     pi = PI
     call RuntimeParameters_get( 'rz_zplane', rz_zplane)
     call RuntimeParameters_get( 'rz_subintNX', rz_subintNX)
     call RuntimeParameters_get( 'rz_subintNY', rz_subintNY)
     call RuntimeParameters_get( 'rz_subintNZ', rz_subintNZ)
     if (NDIM == 2) rz_subintNZ = 1

     !   Parameters for initial particle distribution
#ifdef FLASH_PARTICLES
     call RuntimeParameters_get ("pt_numX", rz_nxp)
     call RuntimeParameters_get( "pt_initialXMin", rz_xpmin)
     call RuntimeParameters_get( "pt_initialXMax", rz_xpmax)
     if (NDIM >= 2) then
        call RuntimeParameters_get ("pt_numY", rz_nyp)
        call RuntimeParameters_get ("pt_initialYMin", rz_ypmin)
        call RuntimeParameters_get( "pt_initialYMax", rz_ypmax)
     else
        rz_nyp = 1
        rz_ypmin = ymin
        rz_ypmax = ymax 
     endif
     if (NDIM == 3) then
        call RuntimeParameters_get ( "pt_numZ", rz_nzp)
        call RuntimeParameters_get ( "pt_initialZMin", rz_zpmin)
        call RuntimeParameters_get ( "pt_initialZMax", rz_zpmax)
     else
        rz_nzp   = 1
        rz_zpmin = zmin
        rz_zpmax = zmax 
     endif
#endif

     !=================================================================================
     !  This section moved from Simulation_initBlock to remove multiple writes
     ! Some checking of the inputs:

     !   Make sure 2 fluids have been defined.

     if (NSPECIES /= 2) then
        call Driver_abortFlash ("Simulation_initBlock:  need two fluids!")
     endif

     !   Make sure 2 or 3 dimensions have been specified.
     if ( NDIM == 1 ) then
        print *, 'Error: NDIM = 1. Problem is not designed for 1d.'
        call Driver_abortFlash('Error: NDIM = 1; require NDIM = 2 or 3.')
     end if

     ! Write a message to stdout describing the problem setup.
     if ( sim_meshMe == MASTER_PE)  then
        if(NDIM == 2) then
           call Logfile_stamp                                & 
                (sim_meshMe, "initializing for 2D shock-cylinder problem", 'Simulation_initBlock')
           write (*,*) "flash:  initializing for 2D shock-cylinder problem"
        elseif(NDIM == 3) then
           call Logfile_stamp                                & 
                (sim_meshMe, "initializing for 3D shock-cylinder problem", 'Simulation_initBlock')
           write (*,*) "flash:  initializing for 3D shock-cylinder problem"
        endif  ! NDIM
        write (*,*)

        if ( ((sim_useRadialFit) .and. (use_rz_sim_data)) .or.        &
             ((sim_useRawData)   .and. (use_rz_sim_data)) .or.        &
             ((sim_useRadialFit) .and. (sim_useRawData)   )     ) then
           write (*,*) "error: Simulation_initBlock: Multiple cylinder initializations defined"
           call Driver_abortFlash ("Simulation_initBlock:  Multiple cylinder initializations defined")
        elseif ( (.NOT. sim_useRadialFit) .and.                 &
             (.NOT. use_rz_sim_data) .and.                 &
             (.NOT. sim_useRawData )      ) then
           write (*,*) "error: Simulation_initBlock: No cylinder initializations defined"
           call Driver_abortFlash ("Simulation_initBlock:  No cylinder initializations defined")
        endif

        if(sim_useRadialFit) then
           call Logfile_stamp                                & 
                (sim_meshMe, "  using radial fit for cylinder initialization", 'Simulation_initBlock')
           write (*,*) "flash:  using radial fit for cylinder initialization"
        endif
        if(use_rz_sim_data) then
           call Logfile_stamp                                & 
                (sim_meshMe, "  using rz simulation data for cylinder initialization", 'Simulation_initBlock')
           write (*,*) "flash:  using rz simulation data for cylinder initialization"
        endif
        if(sim_useRawData) then
           call Logfile_stamp("  using raw experimental image for cylinder initialization", 'Simulation_initBlock')
           write (*,*) "flash:  using raw experimental image for cylinder initialization"
        endif

        write (*,*)
        write (*,*)
     end if

     !=============================================================================================


     !   Read the data files and do preliminary calculations

     ! cc refers to centered in r, centered in z    cell-centered variables
     ! ce refers to centered in r, edge in z        vertical (z) velocity
     ! ec refers to edge in r, centered in z        radial velocity

     ! Data file and domain sizes for rz data.
     call RuntimeParameters_get("nr_c", nr_c )
     call RuntimeParameters_get("nz_c", nz_c )
     if ( (nr_c > nr_c_max) .or. (nz_c > nz_c_max) ) then
        print *, 'Error: sim_rzInitialConditions:  rz data files are too big!'
        call Driver_abortFlash("sim_rzInitialConditions:  rz data files are too big!")
     endif

     nz_e = nz_c - 1
     nr_e = nr_c - 1


     ! Compute coordinate locations for rz data. All locations are on
     !   the interior of the domain. Note there is no data on the boundary.

     dr_rz = rz_rmax/float(nr_c)
     dz_rz = rz_zmax/float(nz_c)
     dri = 1.0/dr_rz
     dzi = 1.0/dz_rz

     ! First centered r is at r=dr_rz/2.
     do ir = 1, nr_c
        r_c(ir) = dr_rz*(0.5 + float(ir-1))
     enddo

     ! First edge r is at r=dr_rz. The centerline is at r=0.
     do ir = 0, nr_e+1
        r_e(ir) = dr_rz*float(ir)
     enddo

     ! First centered z is at z=dz_rz/2.
     do iz = 1, nz_c
        z_c(iz) = dz_rz*(0.5 + float(iz-1))
     enddo

     ! First edge z is at z=dz_rz.
     do iz = 0, nz_e+1
        z_e(iz) = dz_rz*float(iz)
     enddo

     ! Read the files containing sf6 mass fraction, pressure, radial velocity,
     ! and vertical velocity.
     call RuntimeParameters_get("rz_fileSF6_cc", sf6_file )
     call RuntimeParameters_get("rz_filePres_cc", press_file )
     call RuntimeParameters_get("rz_fileRVel_ec", rvel_file )
     call RuntimeParameters_get("rz_fileZVel_ce", zvel_file )

     x_sf6_cc = 0.0
!!$      press_cc = 0.0
     open(unit=70, file=sf6_file, status="old")
!!$      open(unit=71, file=press_file, status="old")
     do iz = 1,nz_c
        read (70,*) (x_sf6_cc(ir,iz), ir=1,nr_c)   ! mass fraction of SF6
!!$        read (71,*) (press_cc(ir,iz), ir=1,nr_c)   ! pressure
     enddo
     close(70)
!!$      close(71)

     rvel_ec = 0.0
     open(unit=72, file=rvel_file, status="old")
     do iz = 1,nz_c
        !       Do not start with 0; data in file starts with ir=1.
        read (72,*) (rvel_ec(ir,iz), ir=1,nr_e)   ! radial velocity
     enddo
     close(72)

     zvel_ce = 0.0
     open(unit=73, file=zvel_file, status="old")
     do iz = 1,nz_e    ! Do not start with 0; data in file starts with iz=1.
        read (73,*) (zvel_ce(ir,iz), ir=1,nr_c)   ! vertical velocity
     enddo
     close(73)

     ! Define local interpolants. Array input sizes differ for cc, ec, and ce.

     !     SF6: centered in r and z.
     call sim_rzInterpolateCC(nr_c,nz_c,dr_rz,dz_rz,               &
          x_sf6_cc,x_sf6_r1,x_sf6_r2,x_sf6_z1,x_sf6_z2)

!!$    !     press: centered in r and z.
!!$      call sim_rzInterpolateCC(nr_c,nz_c,dr_rz,dz_rz,               &
!!$               press_cc,press_r1,press_r2,press_z1,press_z2)

     !     rvel: edge in r, centered in z
     call sim_rzInterpolateEC(nr_c,nz_c,dr_rz,dz_rz,               &
          rvel_ec,rvel_r1,rvel_r2,rvel_z1,rvel_z2)

     !     zvel: centered in r, edge in z
     call sim_rzInterpolateCE(nr_c,nz_c,dr_rz,dz_rz,               &
          zvel_ce,zvel_r1,zvel_r2,zvel_z1,zvel_z2)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  endif  !! use_rz_sim_data


  call RuntimeParameters_get( 'smallx', sim_smallSF6 )

  status = 0
  call Multispecies_getProperty(AIR_SPEC, A, mw_air)
  if (status /= 0)&
       call Driver_abortFlash("Simulation_initBlock:can't find molecular weight of air!")

  call Multispecies_getProperty(SF6_SPEC, A, mw_sf6)
  if (status /= 0)&
       call Driver_abortFlash ("Simulation_initBlock:  can't find molecular weight of sf6!")

  call Multispecies_getProperty(AIR_SPEC, GAMMA,gamma_air)
  if (status /= 0)&
       call Driver_abortFlash ("Simulation_initBlock:  can't find gamma of air!")

  call Multispecies_getProperty(AIR_SPEC, GAMMA, gamma_sf6)
  if (status /= 0)call Driver_abortFlash ("Simulation_initBlock:  can't find gamma of sf6!") 

  !   Initial concentration (mole/cm^3) in ambient air.
  c_initial = rho_amb/mw_air

  !   Z-velocity of sf6.
  !   Assume that only the sf6 has z-momentum. So initialize the z-velocity,
  !     according to the mass fraction of sf6, with vz = vz_sf6 at the 
  !     location of the maximum mass fraction of sf6. At this location
  !     we have maxconc, the mole fraction of sf6. Compute a scaling
  !     factor for the velocity by 
  !        xn_sf6 * vz_fact = vz_sf6
  !     (Later) at each point we will compute xn(isf6), then
  !        vz = xn(isf6)*vz_fact

  !   For 2-d ignore input z-velocity of sf6.
  if (NDIM /= 3)  vz_fact = 0.e0

  if (NDIM == 3) then
     !      compute max sf6 mass fraction from max sf6 mole fraction
     xn_sf6 = maxconc*mw_sf6/  &
          (maxconc*mw_sf6 + (1.e0-maxconc)*mw_air)
     !      at xn_sf6, the z velocity is vz_sf6, so the factor is
     vz_fact = vz_sf6/xn_sf6
  endif

  !   Read the image, if raw data is used for initialization.
  if ( sim_useRawData ) then

     !  Read the image file.
     open(unit=62, file='1cyl_img', status="unknown")

     vmax = 0.0
     do i = 1,sim_rawNumPixelsx
        do j = 1,sim_rawNumPixelsy
           read (62,*) vconc(i,sim_rawNumPixelsy-j+1)
           vmax = max(vmax,vconc(i,sim_rawNumPixelsy-j+1))
        end do
     end do

     close(unit=62,status='keep')

     !       Rescale to from zero to maximum value maxconc.
     do j = 1,sim_rawNumPixelsy
        do i = 1,sim_rawNumPixelsx 
           vconc(i,j) = vconc(i,j)*maxconc/vmax
        end do
     end do

     !       Get dimensions and sizes associated with the image.
     do i = 1,sim_rawNumPixelsx
        sizex(i) = sim_rawMinX + i*sim_rawPixelSize
     end do

     do j = 1,sim_rawNumPixelsy
        sizey(j) = sim_rawMinY + j*sim_rawPixelSize
     end do

     ximgmax = sizex(sim_rawNumPixelsx)
     yimgmax = sizey(sim_rawNumPixelsy)
  end if  ! sim_useRawData

  ! random stuff that Cal must have added after the ShockCyl was moved to FLASH3

  call RuntimeParameters_get('d_pert', d_pert)
  call RuntimeParameters_get('n_pert', n_pert)

  if ( d_pert > 0.e0 ) then

     if ( n_pert < 0 ) then
        n_pert = 0
        write(*,*)          '[Simulation_init] Warning: n_pert reset to 0'
        call Logfile_stamp('[Simulation_init] Warning: n_pert reset to 0')
     end if

     ! randomize the initial velocity perturbations
     call Grid_getNumProcs(NumPEs)

     i = int((n_pert+1) * (sim_meshMe+NumPEs+1))

     call sim_rluxgo (4,i,0,0, (sim_meshMe.eq.MASTER_PE) )

     sdim_i = 1.e0/sqrt(float(NDIM))

  end if
  ! end of blaming Cal.

  species_sf6 = SF6_SPEC - SPECIES_BEGIN + 1
  species_air = AIR_SPEC - SPECIES_BEGIN + 1

  !  Multigamma does not calculate any derivatives, so mask is not needed
  ! vecLen = 1 set in initializations
!!  mask = .false.
  mode = MODE_DENS_PRES




  return

end subroutine Simulation_init
