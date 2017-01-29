!!****if* source/Simulation/SimulationMain/ShockCyl/sim_radialFit
!!
!! NAME
!!
!!  sim_radialFit
!!
!! SYNOPSIS
!!
!!  sim_radialFit(real, intent(in)  :: radius,
!!                real, intent(in)  :: max_sf6_molefrac,
!!                real, intent(out)  :: sf6_molefrac)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   radius : 
!!
!!   max_sf6_molefrac : 
!!
!!   sf6_molefrac : 
!!
!! AUTOGENROBODOC
!!
!!
!!***


!
! Second form of initializing the cylinder region.
!
! Use Tomek's radial fit of the experimental image to set initial
!   SF_6 distribution. Given the distance from the center,
!   return the mole fraction of SF_6.
!
!


subroutine sim_radialFit (radius,max_sf6_molefrac,sf6_molefrac)

!============================================================================

implicit none

real, intent(in)  :: radius
real, intent(in)  :: max_sf6_molefrac
real, intent(out) :: sf6_molefrac

real :: u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12
real :: u13, u14, u15, u16
real :: sim_rawPixelSize, sim_radialNumPixels, sim_radialFitRadius
real :: d, d2, dr, dist, func, maxfunc

integer ::  vuse  ! form of the fitting function
!============================================================================


!   Parameters for fitting function.
    vuse = 3        ! DEV determine by what was set last in original version
    ! v1
    if (vuse == 1) then
        u1 = 274.8318E0
        u2 = 263.7492E0
        u3 = 1.510873E0
        u4 = -2.2534603E-02
        u5 = 1.5714757E-06
        u6 = -3.3182207E-11
        u7 = 171.6132E0
        u8 = 162.8170E0
        u9  = 0.e0
        u10 = 0.e0
        u11 = 0.e0
        u12 = 0.e0
        u13 = 0.e0
        u14 = 0.e0
        u15 = 0.e0
        u16 = 0.e0
     elseif (vuse == 2) then
    ! v2
        u1  = 274.8E0
        u2  = 263.8E0
        u3  = 0.E0
        u4  = 0.E0
        u5  = 0.E0
        u6  = 0.E0
        u7  = 138.9415E0
        u8  = 83.61433E0
        u9  = 11.44282E0
        u10 = 18.68107E0
        u11 = 31.19900E0
        u12 = 39.61253E0
        u13 = 35.13314E0
        u14 = -18.33495E0
        u15 = 87.06977E0
        u16 = 44.60311E0
    elseif (vuse == 3) then
    ! v3
        u1  = 274.8E0
        u2  = 263.8E0
        u3  = 0.E0
        u4  = 0.E0
        u5  = 0.E0
        u6  = 0.E0
        u7  = 144.0725E0
        u8  = 69.45422E0
        u9  = 9.221577E0
        u10 = 20.10299E0
        u11 = 32.47960E0
        u12 = 42.59238E0
        u13 = 32.10067E0
        u14 = -1.559247E0
        u15 = 98.27106E0
        u16 = 15.51837E0
    endif  ! choose fitting function v#
    
    ! the following depends on the fit formula (currently it is common
    !   for all versions)
    maxfunc = u3 + u7 + u9                               &
             + u11*exp(-u12**2/u13**2)                   &
             + u14*exp(-u15**2/u16**2)

    ! Also for all versions (since the same image was fit in the
    !   different versions)

    sim_rawPixelSize = 0.0038   ! cm
    sim_radialNumPixels    = 150

    ! The fit is for 150 pixels radius; if the input radius is more than
    ! 5% bigger than this, just return 0.0 for the mole_frac.
    ! DEV note parameters in Config file not used
    sim_radialFitRadius = 1.05*sim_rawPixelSize*sim_radialNumPixels

!============================================================================


  if(radius <= sim_radialFitRadius) then

    d   = radius/sim_rawPixelSize
    d2  = d*d

    if (vuse .NE. 2) then
        func = u3 + u4*d + u5*d**2 + u6*d**3                 &
                  + u7*exp(-d2/u8**2)                        &
                  + u9*exp(-d2/u10**2)                       &
                  + u11*exp(-(d-u12)**2/u13**2)              &
                  + u14*exp(-(d-u15)**2/u16**2)
    endif ! v1 or v3
    
    ! used only for v2
    if (vuse == 2) then
        dr      = 1.e0 - max(0.e0,min(1.e0,                  &
                  (dist-0.85e0*radius)/(0.85e0*radius) ))
        func    = func * dr
    endif !v2

     sf6_molefrac = max_sf6_molefrac*func/maxfunc

   else  ! radius is outside fit domain

     sf6_molefrac = 0.0

   endif ! radius


   return
 end subroutine sim_radialFit

