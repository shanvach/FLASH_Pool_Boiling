!!****if* source/Simulation/SimulationMain/AdvectMassScalars/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer :: blockId)
!!                       
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up 5 mass scalars in
!!  different shapes to be advected.
!!
!!
!! ARGUMENTS
!!
!!  blockId        The number of the block to initialize
!!
!! PARAMETERS
!!
!!     smallp          smallest pressure allowed
!!     smallx          smallest abundance allowed 
!!     rhoin           Density inside the pulse
!!     rhoout          Density outside the pulse
!!     pressure        Pressure
!!     velocity        Fluid velocity
!!     posn            Position of the pulse center at x-axis (y=z=0)
!!     width           Width of the pulse along x-axis
!!     phase           Phase shift between density and tracers
!!     xangle          Angle made by ndiaphragm normal w/x-axis (deg)
!!     yangle          Angle made by diaphragm normal w/y-axis (deg)
!!     pulse_fctn      Which pulse shape function to use
!!
!!
!!***

subroutine Simulation_initBlock(blockId)

  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getBlkPhysicalSize, Grid_getBlkCenterCoords, Grid_getDeltas, &
    Grid_putPointData
  
  implicit none

#include "constants.h"
#include "Flash.h"

  integer,intent(IN) ::  blockId
  

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
!!$  integer, parameter :: q = max(nxb+2*nguard, & 
!!$       &                              nyb+2*nguard*k2d, & 
!!$       &                              nzb+2*nguard*k3d)

  integer :: i, j, k, n
  integer :: imax, jmax, kmax

  integer :: ii, jj, kk, Nint
  real    :: Nintinv1, NintVol_i

  real :: sum_rho, sum_ms1, sum_ms2, sum_ms3, sum_ms4, sum_ms5

  real :: delx, xx, dely, yy, delz, zz
  real :: lposn0, lposn
  real :: xpos, ypos, zpos

  real :: xbmin, xbmax, ybmin, ybmax, zbmin, zbmax

  integer :: blkLimits(HIGH, MDIM), blkLimitsGC(HIGH, MDIM)
  integer :: axis(MDIM)

  real :: rho_zone, velx_zone, vely_zone, velz_zone, pres_zone, & 
       &     ener_zone, ekin_zone

  real :: phase

  real :: del(MDIM)
  real :: blockSize(MDIM), blockCenter(MDIM)
  real :: wfac
  real :: PulseShape
  real :: ms1, ms2, ms3, ms4, ms5

  !==============================================================================

  call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC)

  imax = blkLimitsGC(HIGH, IAXIS)    ! Maximum values of the cell index ranges
  jmax = blkLimitsGC(HIGH, JAXIS)
  kmax = blkLimitsGC(HIGH, KAXIS)

  call Grid_getBlkPhysicalSize(blockId, blockSize)
  call Grid_getBlkCenterCoords(blockId, blockCenter)

  ! compute the extrema of the current blocks
  xbmax = blockCenter(IAXIS) + 0.5e0*blockSize(IAXIS)
  xbmin = blockCenter(IAXIS) - 0.5e0*blockSize(IAXIS)

  ybmax = blockCenter(JAXIS) + 0.5e0*blockSize(JAXIS)
  ybmin = blockCenter(JAXIS) - 0.5e0*blockSize(JAXIS)

  zbmax = blockCenter(KAXIS) + 0.5e0*blockSize(KAXIS)
  zbmin = blockCenter(KAXIS) - 0.5e0*blockSize(KAXIS)

  ! compute the cell size
  call Grid_getDeltas(blockId, del)
  
  delx = del(IAXIS)
  dely = del(JAXIS)
  delz = del(KAXIS)

  Nint = 5
  Nintinv1  = 1.e0/(float(Nint) + 1.e0)
  NintVol_i = 1.e0/(float(Nint)**NDIM)

  !=============================================================================

  ! Loop over cells in the block.  For each, compute the physical position of 
  ! its left and right edge and its center as well as its physical width.  Then 
  ! decide which side of the initial discontinuity it is on and initialize the 
  ! hydro variables appropriately.

  ! Dongwook - phase has not been initialized.
  phase = 0.

  do k = 1, kmax
     do j = 1, jmax
        do i = 1, imax

           ! For each zone, find the average density by computing the pulse shape at a 
           ! number of points inside the zone.

           sum_rho = 0.e0
           sum_ms1 = 0.e0
           sum_ms2 = 0.e0
           sum_ms3 = 0.e0
           sum_ms4 = 0.e0
           sum_ms5 = 0.e0

           do kk = 1, 1+(Nint-1)*K3D

              ! compute the subzone z coordinate
              zz = zbmin + delz*(float(k-NGUARD-1)+kk*Nintinv1)
              lposn0 = sim_posn - zz*sim_zcos/sim_xcos
              zpos = zz
!print*,zz,lposn0

              do jj = 1, 1+(Nint-1)*K2D

                 ! compute the subzone y coordinate
                 yy = ybmin + dely*(float(j-NGUARD-1)+jj*Nintinv1)
                 lposn = lposn0 - yy*sim_ycos/sim_xcos
                 ypos = yy * sim_ycos

                 do ii = 1, Nint

                    ! compute the subzone x coordinate
                    xx  = xbmin + delx*(float(i-NGUARD-1)+ii*Nintinv1)
                    xpos = (xx - lposn) * sim_xcos / sim_width

                    ! compute the weighting fraction for the pulse shape, and use this to
                    ! weight the current subzone's contribution to the current zone
#if NDIM == 2
                    if ( sim_pulse_fctn == 3 .or. sim_pulse_fctn == 4 .or. sim_pulse_fctn == 5 ) then
                       call sim_pulseShape(xx,   yy  , sim_pulse_fctn, 0.e0, wfac)
                       sum_rho = sum_rho + sim_rhoin*wfac + sim_rhoin + sim_rhoout
                    else
                       call sim_pulseShape(xpos, ypos, sim_pulse_fctn, 0.e0, wfac)
                       sum_rho = sum_rho + sim_rhoin*wfac + sim_rhoout*(1.e0-wfac) 
                    end if

                    call sim_pulseShape(xpos, ypos, sim_pulse_fctn_ms1, phase, wfac)

                    sum_ms1 = sum_ms1 + sim_msin*wfac + sim_msout*(1.e0-wfac) 

                    call sim_pulseShape(xx, yy, sim_pulse_fctn_ms2, phase, wfac)

                    sum_ms2 = sum_ms2 + sim_msin*wfac + sim_msout*(1.e0-wfac) 

                    call sim_pulseShape(xx, yy, sim_pulse_fctn_ms3, phase, wfac)

                    sum_ms3 = sum_ms3 + sim_msin*wfac + sim_msout*(1.e0-wfac) 

                    call sim_pulseShape(xx, yy, sim_pulse_fctn_ms4, phase, wfac)

                    sum_ms4 = sum_ms4 + sim_msin*wfac + sim_msin + sim_msout

                    call sim_pulseShape(xx, yy, sim_pulse_fctn_ms5, phase, wfac)

                    sum_ms5 = sum_ms5 + sim_msin*wfac

#elif NDIM == 3

                    if ( sim_pulse_fctn == 3 .or. sim_pulse_fctn == 4 .or. sim_pulse_fctn == 5 ) then
                       call sim_pulseShape3D(xx,   yy  , zz,  sim_pulse_fctn, 0.e0, wfac)
                       sum_rho = sum_rho + sim_rhoin*wfac + sim_rhoin + sim_rhoout
                    else
                       call sim_pulseShape3D(xpos, ypos, zpos, sim_pulse_fctn, 0.e0, wfac)
                       sum_rho = sum_rho + sim_rhoin*wfac + sim_rhoout*(1.e0-wfac) 
                    end if

                    call sim_pulseShape3D(xpos, ypos, zpos, sim_pulse_fctn_ms1, phase, wfac)

                    sum_ms1 = sum_ms1 + sim_msin*wfac + sim_msout*(1.e0-wfac) 

                    call sim_pulseShape3D(xx, yy, zz, sim_pulse_fctn_ms2, phase, wfac)

                    sum_ms2 = sum_ms2 + sim_msin*wfac + sim_msout*(1.e0-wfac) 

                    call sim_pulseShape3D(xx, yy, zz, sim_pulse_fctn_ms3, phase, wfac)

                    sum_ms3 = sum_ms3 + sim_msin*wfac + sim_msout*(1.e0-wfac) 

                    call sim_pulseShape3D(xx, yy, zz,sim_pulse_fctn_ms4, phase, wfac)

                    sum_ms4 = sum_ms4 + sim_msin*wfac + sim_msin + sim_msout

                    call sim_pulseShape3D(xx, yy,zz,  sim_pulse_fctn_ms5, phase, wfac)

                    sum_ms5 = sum_ms5 + sim_msin*wfac
#endif

                 end do

              end do

           end do

           ! Initialize the hydro quantities.
           rho_zone = sum_rho * NintVol_i 
!!$           rho_zone = 1.

           ms1      = sum_ms1 * NintVol_i 
           ms2      = sum_ms2 * NintVol_i 
           ms3      = sum_ms3 * NintVol_i 
           ms4      = sum_ms4 * NintVol_i 
           ms5      = sum_ms5 * NintVol_i 

           pres_zone = sim_pressure
!!$           pres_zone = 1. 


           velx_zone = sim_velocity * sim_xcos
           vely_zone = sim_velocity * sim_ycos
           velz_zone = sim_velocity * sim_zcos

           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k

#if NSPECIES > 0
           ! initialize the nuclear abundances
           Call Grid_putPointData(blockId, CENTER, SPECIES_BEGIN, EXTERIOR, axis, 1.0e0-(NSPECIES-1)*sim_smallX)

           do n = SPECIES_BEGIN+1,SPECIES_END
              call Grid_putPointData(blockId, CENTER, n, EXTERIOR, axis, sim_smallX)
           end do
#endif

!           if ( ims1 > 0 ) call dBasePutData(ims1, iPt, i, j, k, block_no, ms1, 1)
           if ( sim_ims1 > 0 ) call Grid_putPointData(blockId, CENTER, sim_ims1, EXTERIOR, axis, ms1)
           if ( sim_ims2 > 0 ) call Grid_putPointData(blockId, CENTER, sim_ims2, EXTERIOR, axis, ms2)
           if ( sim_ims3 > 0 ) call Grid_putPointData(blockId, CENTER, sim_ims3, EXTERIOR, axis, ms3)
           if ( sim_ims4 > 0 ) call Grid_putPointData(blockId, CENTER, sim_ims4, EXTERIOR, axis, ms4)
           if ( sim_ims5 > 0 ) call Grid_putPointData(blockId, CENTER, sim_ims5, EXTERIOR, axis, ms5)

           ! compute the gas energy and set the gamma-values needed for the equation of 
           ! state.

           ekin_zone = 0.5e0*( velx_zone**2 & 
                              +vely_zone**2 & 
                              +velz_zone**2)

           ener_zone = pres_zone / (sim_gamma-1.e0)
           ener_zone = ener_zone / rho_zone
           ener_zone = ener_zone + ekin_zone
           ener_zone = max(ener_zone, sim_smallp)

           ! update the variables in the current zone via the database put method

           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho_zone)
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, pres_zone)
           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, ener_zone)

           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, sim_gamma)
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, sim_gamma)

           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velx_zone)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, vely_zone)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, velz_zone)

        end do
     end do
  end do

  return
end subroutine Simulation_initBlock
