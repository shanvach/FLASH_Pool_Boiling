!!****if* source/Simulation/SimulationMain/ShockCyl/sim_rzInitialConditions
!!
!! NAME
!!
!!  sim_rzInitialConditions
!!
!! SYNOPSIS
!!
!!  sim_rzInitialConditions(real, intent(in)  :: radius,
!!                          real, intent(in)  :: height,
!!                          real, intent(out)  :: sf6_massfrac,
!!                          real, intent(out)  :: rvel,
!!                          real, intent(out)  :: zvel)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   radius : 
!!
!!   height : 
!!
!!   sf6_massfrac : 
!!
!!   rvel : 
!!
!!   zvel : 
!!
!! AUTOGENROBODOC
!!
!!
!!***


!
! Return interpolated values for sf6_massfrac, press, rvel, and zvel
!   given radius and height.
!
! In the first call block, data files are read which contain data on
! an equispaced r-z grid
!   sf6_massfrac    centered in r, centered in z
!   press           centered in r, centered in z
!   rvel            edge     in r, centered in z
!   zvel            centered in r, edge in z
! Interpolants for each variable are defined, and the interpolants
!   have slightly different domains depending on how the data is centered.
!
! After the first call, the variable values are computed from the
!   interpolants for the given radius and height.
!
! There is no limiting applied and no checking of the datafiles.
! If limiting (of e.g. press > 0.0) is desired, this should be handled
!   in the calling program.
!

!  Drop pressure - since is is essentially constant, just use p_ambient
!    in the calling program.
!


subroutine sim_rzInitialConditions(radius, height, sf6_massfrac, rvel, zvel)


!  use runtime_parameters, ONLY: get_parm_from_context,GLOBAL_PARM_CONTEXT
  use Simulation_data
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none 

#include "constants.h"
!  include 'sim_rzDatafileSizes.fh'


  real, intent(in)  :: radius, height
  real, intent(out) :: sf6_massfrac, rvel, zvel

!!$  real, save    :: press_cc(  nr_c_max,   nz_c_max),  &
!!$                   press_r1(  nr_c_max,   nz_c_max),  &
!!$                   press_r2(  nr_c_max,   nz_c_max),  &
!!$                   press_z1(  nr_c_max,   nz_c_max),  &
!!$                   press_z2(  nr_c_max,   nz_c_max)





!============================================================================

! Some input checking.

  if(radius > rz_rmax) then
    print *,  'Error: sim_rzInitialConditions: Interpolants are not defined at this radius,',radius
    call Driver_abortFlash("sim_rzInitialConditions:: Interpolants are not defined at this radius.")
  elseif(height < 0.0) then 
    print *,  'Error: sim_rzInitialConditions:: Interpolants are not defined at this height,',height
    call Driver_abortFlash("sim_rzInitialConditions:: Interpolants are not defined at this height.")
  elseif(height > rz_zmax) then 
    print *,  'Error: sim_rzInitialConditions:: Interpolants are not defined at this height,',height
    call Driver_abortFlash("sim_rzInitialConditions:: Interpolants are not defined at this height.")
  endif

! Given the radius, compute the sf6 mass fraction, pressure, radial
!   velocity and vertical velocity using the interpolants defined
!   in the first_call block.


! Find the centered and edge indices nearest to the input radius and height.
! Recall r(ir_c) = dr/2 + (ir_c-1)*dr
!        r(ir_e) = (ir_e)*dr
  ir_c = 1 + int(radius/dr_rz)
  ir_e = ir_c
  if ( (mod(radius,dr_rz)/dr_rz) < 0.5 ) ir_e = ir_e - 1

  iz_c = 1 + int(height/dz_rz)
  iz_e = iz_c
  if ( (mod(height,dz_rz)/dz_rz) < 0.5 ) iz_e = iz_e - 1

  sf6_massfrac =   x_sf6_cc(ir_c,iz_c)                            &
                 + x_sf6_r1(ir_c,iz_c)*( radius - r_c(ir_c) )     &
                 + x_sf6_r2(ir_c,iz_c)*( radius - r_c(ir_c) )**2  &
                 + x_sf6_z1(ir_c,iz_c)*( height - z_c(iz_c) )     &
                 + x_sf6_z2(ir_c,iz_c)*( height - z_c(iz_c) )**2

!!$  press        =   press_cc(ir_c,iz_c)                            &
!!$                 + press_r1(ir_c,iz_c)*( radius - r_c(ir_c) )     &
!!$                 + press_r2(ir_c,iz_c)*( radius - r_c(ir_c) )**2  &
!!$                 + press_z1(ir_c,iz_c)*( height - z_c(iz_c) )     &
!!$                 + press_z2(ir_c,iz_c)*( height - z_c(iz_c) )**2

!   No  bounding  0.0 < sf6 < 1.0
!                 0.0 < press
!   The calling program should limit appropriately 

   rvel        =   rvel_ec(ir_e,iz_c)                            &
                 + rvel_r1(ir_e,iz_c)*( radius - r_e(ir_e) )     &
                 + rvel_r2(ir_e,iz_c)*( radius - r_e(ir_e) )**2  &
                 + rvel_z1(ir_e,iz_c)*( height - z_c(iz_c) )     &
                 + rvel_z2(ir_e,iz_c)*( height - z_c(iz_c) )**2

   zvel        =   zvel_ce(ir_c,iz_e)                            &
                 + zvel_r1(ir_c,iz_e)*( radius - r_c(ir_c) )     &
                 + zvel_r2(ir_c,iz_e)*( radius - r_c(ir_c) )**2  &
                 + zvel_z1(ir_c,iz_e)*( height - z_e(iz_e) )     &
                 + zvel_z2(ir_c,iz_e)*( height - z_e(iz_e) )**2

  return
end subroutine sim_rzInitialConditions

subroutine sim_rzInterpolateCC(nr_c,nz_c,dr,dz,cc_data,cc_r1,cc_r2,cc_z1,cc_z2)

!
! Define 2d quadratic interpolants on a 5-point stencil.
!
! Assumed grid is in r, z with equal spacing dr, dz.
!    0 < r < rmax,    nr_c intervals
!    0 < z < zmax,    nz_c intervals
! Data is at cell centers, with
!   ir=1     at r=dr/2
!   ir=nr_c  at r=dr/2 + (nr_c-1)*dr
!   iz=1     at z=dz/2
!   iz=nz_c  at z=dz/2 + (nz_c-1)*dz
!
! Boundary conditions, assuming thermodynamic variable q:
!   r=0      dq/dr=0   => q(ir=0) = q(ir=1)
!   r=rmax   one-sided linear interpolation on interior, d2q/dr2 = 0
!   z=0      dq/dz=0   => q(iz=0) = q(iz=1)
!   z=zmax   dq/dz=0   => q(iz=nz_c+1) = q(iz=nz_c)
!
! Then for a point (r,z)  nearest to (r(ir), z(iz))
!    q(r,z) = q(ir,iz) + (r-r(ir))*cc_r1(ir,iz) + (r-r(ir))^2*cc_r2(ir,iz)
!                      + (z-z(iz))*cc_z1(ir,iz) + (z-z(iz))^2*cc_z2(ir,iz)
!
!
  implicit none

  include 'sim_rzDatafileSizes.fh'

  integer, intent(in)    :: nr_c, nz_c
  real,    intent(in)    :: dr, dz
  real,    intent(in)    :: cc_data(nr_c_max, nz_c_max)
  real,    intent(inout) :: cc_r1(nr_c_max, nz_c_max),  &
                            cc_r2(nr_c_max, nz_c_max),  &
                            cc_z1(nr_c_max, nz_c_max),  &
                            cc_z2(nr_c_max, nz_c_max)
  integer  :: ir, iz
  real     :: dri, dzi

! Define local interpolants.

  dri = 1.0e0/dr
  dzi = 1.0e0/dz

  do iz = 2, nz_c-1
!   Near centerline (ir=1)
    cc_r1(1 ,iz) = 0.5*dzi*(  cc_data(2  ,iz )              &
                            - cc_data(1  ,iz ) )
    cc_z1(1 ,iz) = 0.5*dzi*(  cc_data(1  ,iz+1)               &
                            - cc_data(1  ,iz-1)) 
    cc_r2(1 ,iz) = 0.5*dri**2*(      cc_data(2   ,iz  )     &
                                   - cc_data(1   ,iz  ) ) 
    cc_z2(1 ,iz) = 0.5*dzi**2*(      cc_data(1   ,iz+1)       &
                               - 2.0*cc_data(1   ,iz  )       &
                                   + cc_data(1   ,iz-1) ) 
!   Interior:
    do ir = 2, nr_c-1
        cc_r1(ir,iz) = 0.5*dri*(  cc_data(ir+1,  iz)          &
                                - cc_data(ir-1,  iz) ) 
        cc_z1(ir,iz) = 0.5*dzi*(  cc_data(ir  ,iz+1)          &
                                - cc_data(ir  ,iz-1)) 
        cc_r2(ir,iz) = 0.5*dri**2*(      cc_data(ir+1,iz  )   &
                                   - 2.0*cc_data(ir  ,iz  )   &
                                       + cc_data(ir-1,iz  ) ) 
        cc_z2(ir,iz) = 0.5*dzi**2*(      cc_data(ir  ,iz+1)   &
                                   - 2.0*cc_data(ir  ,iz  )   &
                                       + cc_data(ir  ,iz-1) ) 
    enddo
!   Far field (ir=nr_c). Assume linear fit between nr_c-1 and nr_c
    cc_r1(nr_c,iz) = dri*(  cc_data(nr_c  ,  iz)            &
                          - cc_data(nr_c-1,  iz) ) 
    cc_z1(nr_c,iz) = 0.5*dzi*(  cc_data(nr_c  ,iz+1)          &
                              - cc_data(nr_c  ,iz-1) ) 
    cc_r2(nr_c,iz) = 0.0
    cc_z2(nr_c,iz) = 0.5*dzi**2*(      cc_data(nr_c  ,iz+1)   &
                                 - 2.0*cc_data(nr_c  ,iz  )   &
                                     + cc_data(nr_c  ,iz-1) )
  enddo

! Bottom wall (iz=1): Reflection in z.
  do ir = 2, nr_c-1
    cc_r1(ir,1 ) = 0.5*dri*(  cc_data(ir+1,  1 )              &
                            - cc_data(ir-1,  1 ) ) 
    cc_z1(ir,1 ) = 0.5*dzi*(  cc_data(ir  ,  2 )            &
                            - cc_data(ir  ,  1 ) )
    cc_r2(ir,1 ) = 0.5*dri**2*(      cc_data(ir+1,1   )       &
                               - 2.0*cc_data(ir  ,1   )       &
                                   + cc_data(ir-1,1   ) ) 
    cc_z2(ir,1 ) = 0.5*dzi**2*(  cc_data(ir  ,2   )         &
                               - cc_data(ir  ,1   ) ) 
  enddo

! Top wall (iz=nz_c): Reflection in z.
  do ir = 2, nr_c-1
    cc_r1(ir,nz_c ) = 0.5*dri*(  cc_data(ir+1,nz_c )          &
                               - cc_data(ir-1,nz_c ) ) 
    cc_z1(ir,nz_c ) = 0.5*dzi*(  cc_data(ir  ,nz_c  )      &
                               - cc_data(ir  ,nz_c-1) )
    cc_r2(ir,nz_c ) = 0.5*dri**2*(      cc_data(ir+1,nz_c)    &
                                  - 2.0*cc_data(ir  ,nz_c)    &
                                      + cc_data(ir-1,nz_c) ) 
    cc_z2(ir,nz_c ) = 0.5*dzi**2*(- cc_data(ir  ,nz_c  )   &
                                  + cc_data(ir  ,nz_c-1) ) 
  enddo

! Corners.
  cc_r1(1,1) = 0.5*dri*(  cc_data(2  , 1)    &
                        - cc_data(1  , 1) )
  cc_z1(1,1) = 0.5*dzi*(  cc_data(1  , 2)    &
                        - cc_data(1  , 1) )
  cc_r2(1,1) = 0.5*dri**2*(  cc_data(2  ,1  )   &
                           - cc_data(1  ,1  ) ) 
  cc_z2(1,1) = 0.5*dzi**2*(  cc_data(1  ,2  )   &
                           - cc_data(1  ,1  ) ) 


  cc_r1(1,nz_c) = 0.5*dri*(  cc_data(2  ,nz_c)      &
                           - cc_data(1  ,nz_c) )
  cc_z1(1,nz_c) = 0.5*dzi*(  cc_data(1  ,nz_c  )    &
                           - cc_data(1  ,nz_c-1) )
  cc_r2(1,nz_c) = 0.5*dri**2*(  cc_data(2  ,nz_c)     &
                              - cc_data(1  ,nz_c) ) 
  cc_z2(1,nz_c) = 0.5*dzi**2*(- cc_data(1  ,nz_c  )   &
                              + cc_data(1  ,nz_c-1) ) 


  cc_r1(nr_c,1) = dri*(  cc_data(nr_c  ,  1)     &
                       - cc_data(nr_c-1,  1) )
  cc_z1(nr_c,1) = 0.5*dzi*(  cc_data(nr_c, 2)    &
                           - cc_data(nr_c, 1) )
  cc_r2(nr_c,1) = 0.0
  cc_z2(nr_c,1) = 0.5*dzi**2*(  cc_data(nr_c,2  )   &
                              - cc_data(nr_c,1  ) ) 



  cc_r1(nr_c,nz_c) = dri*(  cc_data(nr_c  ,nz_c)        &
                          - cc_data(nr_c-1,nz_c) )
  cc_z1(nr_c,nz_c) = 0.5*dzi*(  cc_data(nr_c,nz_c  )    &
                              - cc_data(nr_c,nz_c-1) )
  cc_r2(nr_c,nz_c) = 0.0
  cc_z2(nr_c,nz_c) = 0.5*dzi**2*(- cc_data(nr_c,nz_c  )   &
                                 + cc_data(nr_c,nz_c-1) ) 

  return

end subroutine sim_rzInterpolateCC

subroutine sim_rzInterpolateCE(nr_c,nz_c,dr,dz,ce_data,ce_r1,ce_r2,ce_z1,ce_z2)

!
! Define 2d quadratic interpolants on a 5-point stencil.
!
! Data is assumed to be vertical velocity.
!
! Assumed grid is in r, z with equal spacing dr, dz.
!    0 < r < rmax,    nr_c intervals
!    0 < z < zmax,    nz_c intervals
! Data is at cell centers in r, edges in z, with
!   ir=1     at r=dr/2
!   ir=nr_c  at r=dr/2 + (nr_c-1)*dr
!
!   iz=0     at z=0
!   iz=1     at z=dz
!   iz=nz_c  at z=zmax= (nz_c)*dz
!
! Boundary conditions, assuming vertical velocity q:
!   Centerline:
!   r=0      dq/dr=0   => q(ir=0) = q(ir=1)
!   r=rmax   one-sided linear interpolation on interior, d2q/dr2 = 0
!   Slip walls, upper and lower. No flow through.
!   z=0      q=0 and dq/dz=0 
!   z=zmax   q=0 and dq/dz=0
!
! Then for a point (r,z)  nearest to (r(ir), z(iz))
!    q(r,z) = q(ir,iz) + (r-r(ir))*ce_r1(ir,iz) + (r-r(ir))^2*ce_r2(ir,iz)
!                      + (z-z(iz))*ce_z1(ir,iz) + (z-z(iz))^2*ce_z2(ir,iz)
!
!
  implicit none
  include 'sim_rzDatafileSizes.fh'
  integer, intent(in)    :: nr_c, nz_c
  real,    intent(in)    :: dr, dz
  real,    intent(in)    :: ce_data(nr_c_max, 0:nz_c_max)
  real,    intent(inout) :: ce_r1(nr_c_max, 0:nz_c_max),  &
                            ce_r2(nr_c_max, 0:nz_c_max),  &
                            ce_z1(nr_c_max, 0:nz_c_max),  &
                            ce_z2(nr_c_max, 0:nz_c_max)
  integer  :: ir, iz
  real     :: dri, dzi

! Define local interpolants.

  dri = 1.0e0/dr
  dzi = 1.0e0/dz

  do iz = 1, nz_c-1
!   Near centerline (ir=1): reflection in r.
    ce_r1(1 ,iz) = 0.5*dzi*(  ce_data(2  ,iz )              &
                            - ce_data(1  ,iz ) )
    ce_z1(1 ,iz) = 0.5*dzi*(  ce_data(1  ,iz+1)               &
                            - ce_data(1  ,iz-1)) 
    ce_r2(1 ,iz) = 0.5*dri**2*(      ce_data(2   ,iz  )     &
                                   - ce_data(1   ,iz  ) ) 
    ce_z2(1 ,iz) = 0.5*dzi**2*(      ce_data(1   ,iz+1)       &
                               - 2.0*ce_data(1   ,iz  )       &
                                   + ce_data(1   ,iz-1) ) 
!   Interior:
    do ir = 2, nr_c-1
        ce_r1(ir,iz) = 0.5*dri*(  ce_data(ir+1,  iz)          &
                                - ce_data(ir-1,  iz) ) 
        ce_z1(ir,iz) = 0.5*dzi*(  ce_data(ir  ,iz+1)          &
                                - ce_data(ir  ,iz-1)) 
        ce_r2(ir,iz) = 0.5*dri**2*(      ce_data(ir+1,iz  )   &
                                   - 2.0*ce_data(ir  ,iz  )   &
                                       + ce_data(ir-1,iz  ) ) 
        ce_z2(ir,iz) = 0.5*dzi**2*(      ce_data(ir  ,iz+1)   &
                                   - 2.0*ce_data(ir  ,iz  )   &
                                       + ce_data(ir  ,iz-1) ) 
    enddo
!   Far field (ir=nr_c). Assume linear fit between nr_c-1 and nr_c
    ce_r1(nr_c,iz) = dri*(  ce_data(nr_c  ,  iz)            &
                          - ce_data(nr_c-1,  iz) ) 
    ce_z1(nr_c,iz) = 0.5*dzi*(  ce_data(nr_c  ,iz+1)          &
                              - ce_data(nr_c  ,iz-1) ) 
    ce_r2(nr_c,iz) = 0.0
    ce_z2(nr_c,iz) = 0.5*dzi**2*(      ce_data(nr_c  ,iz+1)   &
                                 - 2.0*ce_data(nr_c  ,iz  )   &
                                     + ce_data(nr_c  ,iz-1) )
  enddo

! Bottom wall (iz=1): No flow through in z.
  do ir = 2, nr_c-1
    ce_r1(ir,0 ) = 0.5*dri*(  ce_data(ir+1,  0 )              &
                            - ce_data(ir-1,  0 ) ) 
    ce_z1(ir,0 ) =     dzi*(  ce_data(ir  ,  1 ) )
    ce_r2(ir,0 ) = 0.5*dri**2*(      ce_data(ir+1,0   )       &
                               - 2.0*ce_data(ir  ,0   )       &
                                   + ce_data(ir-1,0   ) ) 
    ce_z2(ir,0 ) = 0.0 
  enddo

! Top wall (iz=nz_c): Reflection in z.
  do ir = 2, nr_c-1
    ce_r1(ir,nz_c ) = 0.5*dri*(  ce_data(ir+1,nz_c )          &
                               - ce_data(ir-1,nz_c ) ) 
    ce_z1(ir,nz_c ) =     dzi*(- ce_data(ir  ,nz_c-1) ) 
    ce_r2(ir,nz_c ) = 0.5*dri**2*(      ce_data(ir+1,nz_c)    &
                                  - 2.0*ce_data(ir  ,nz_c)    &
                                      + ce_data(ir-1,nz_c) ) 
    ce_z2(ir,nz_c ) = 0.0
  enddo

! Corners.
  ce_r1(1,0) = 0.5*dri*(  ce_data(2  , 0)    &
                        - ce_data(1  , 0) )
  ce_z1(1,0) =     dzi*(  ce_data(1  , 1  ) )
  ce_r2(1,0) = 0.5*dri**2*(  ce_data(2  ,0  )   &
                           - ce_data(1  ,0  ) ) 
  ce_z2(1,0) = 0.0


  ce_r1(1,nz_c) = 0.5*dri*(  ce_data(2  ,nz_c)      &
                           - ce_data(1  ,nz_c) )
  ce_z1(1,nz_c) =     dzi*(- ce_data(1  ,nz_c-1) )
  ce_r2(1,nz_c) = 0.5*dri**2*(  ce_data(2  ,nz_c)     &
                              - ce_data(1  ,nz_c) ) 
  ce_z2(1,nz_c) = 0.0


  ce_r1(nr_c,0) = dri*(  ce_data(nr_c  ,  0)     &
                       - ce_data(nr_c-1,  0) )
  ce_z1(nr_c,0) = dzi*(  ce_data(nr_c  ,1  ) )
  ce_r2(nr_c,0) = 0.0
  ce_z2(nr_c,0) = 0.0



  ce_r1(nr_c,nz_c) = dri*(  ce_data(nr_c  ,nz_c)        &
                          - ce_data(nr_c-1,nz_c) )
  ce_z1(nr_c,nz_c) = dzi*(- ce_data(nr_c  ,nz_c-1) )
  ce_r2(nr_c,nz_c) = 0.0
  ce_z2(nr_c,nz_c) = 0.0

  return

end subroutine sim_rzInterpolateCE

subroutine sim_rzInterpolateEC(nr_c,nz_c,dr,dz,ec_data,ec_r1,ec_r2,ec_z1,ec_z2)

!
! Define 2d quadratic interpolants on a 5-point stencil.
!
! Data is assumed to be radial velocity.
!
! Assumed grid is in r, z with equal spacing dr, dz.
!    0 < r < rmax,    nr_c intervals
!    0 < z < zmax,    nz_c intervals
! Data is at cell edges in r, centers in z, with
!   ir=0     at r=0
!   ir=1     at r=dr
!   ir=nr_c  at r=rmax= (nr_c)*dr
!
!   iz=1     at z=dz/2
!   iz=nz_c  at z=dz/2 + (nz_c-1)*dz
!
! Boundary conditions, assuming radial velocity q:
!   Centerline.
!   r=0      q=0 and dq/dr=0
!   Far field. 
!   r=rmax   one-sided linear interpolation on interior, d2q/dr2 = 0
!   Slip walls upper and lower.
!   z=0      dq/dz=0   => q(iz=0) = q(iz=1)
!   z=zmax   dq/dz=0   => q(iz=nz_c+1) = q(iz=nz_c)
!
! Then for a point (r,z)  nearest to (r(ir), z(iz))
!    q(r,z) = q(ir,iz) + (r-r(ir))*ec_r1(ir,iz) + (r-r(ir))^2*ec_r2(ir,iz)
!                      + (z-z(iz))*ec_z1(ir,iz) + (z-z(iz))^2*ec_z2(ir,iz)
!
!
  implicit none
  include 'sim_rzDatafileSizes.fh'
  integer, intent(in)    :: nr_c, nz_c
  real,    intent(in)    :: dr, dz
  real,    intent(in)    :: ec_data(0:nr_c_max, nz_c_max)
  real,    intent(inout) :: ec_r1(0:nr_c_max, nz_c_max),   &
                            ec_r2(0:nr_c_max, nz_c_max),   &
                            ec_z1(0:nr_c_max, nz_c_max),   &
                            ec_z2(0:nr_c_max, nz_c_max)
  integer  :: ir, iz
  real     :: dri, dzi

! Define local interpolants.

  dri = 1.0e0/dr
  dzi = 1.0e0/dz

  do iz = 2, nz_c-1
!   Near centerline (ir=1)
    ec_r1(0 ,iz) =     dri*(  ec_data(1,  iz  ) )
    ec_z1(0 ,iz) = 0.5*dzi*(  ec_data(0  ,iz+1)               &
                            - ec_data(0  ,iz-1)) 
    ec_r2(0 ,iz) = 0.0
    ec_z2(0 ,iz) = 0.5*dzi**2*(      ec_data(0   ,iz+1)       &
                               - 2.0*ec_data(0   ,iz  )       &
                                   + ec_data(0   ,iz-1) ) 
!   Interior:
    do ir = 1, nr_c-1
        ec_r1(ir,iz) = 0.5*dri*(  ec_data(ir+1,  iz)          &
                                - ec_data(ir-1,  iz) ) 
        ec_z1(ir,iz) = 0.5*dzi*(  ec_data(ir  ,iz+1)          &
                                - ec_data(ir  ,iz-1)) 
        ec_r2(ir,iz) = 0.5*dri**2*(      ec_data(ir+1,iz  )   &
                                   - 2.0*ec_data(ir  ,iz  )   &
                                       + ec_data(ir-1,iz  ) ) 
        ec_z2(ir,iz) = 0.5*dzi**2*(      ec_data(ir  ,iz+1)   &
                                   - 2.0*ec_data(ir  ,iz  )   &
                                       + ec_data(ir  ,iz-1) ) 
    enddo
!   Far field (ir=nr_c). Assume linear fit between nr_c-1 and nr_c
    ec_r1(nr_c,iz) = dri*(- ec_data(nr_c-1,  iz) ) 
    ec_z1(nr_c,iz) = 0.5*dzi*(  ec_data(nr_c  ,iz+1)          &
                              - ec_data(nr_c  ,iz-1) ) 
    ec_r2(nr_c,iz) = 0.0
    ec_z2(nr_c,iz) = 0.5*dzi**2*(      ec_data(nr_c  ,iz+1)   &
                                 - 2.0*ec_data(nr_c  ,iz  )   &
                                     + ec_data(nr_c  ,iz-1) )
  enddo

! Bottom wall (iz=1): Reflection in z.
  do ir = 1, nr_c-1
    ec_r1(ir,1 ) = 0.5*dri*(  ec_data(ir+1,  1 )              &
                            - ec_data(ir-1,  1 ) ) 
    ec_z1(ir,1 ) = 0.5*dzi*(  ec_data(ir  ,  2 )            &
                            - ec_data(ir  ,  1 ) )
    ec_r2(ir,1 ) = 0.5*dri**2*(      ec_data(ir+1,1   )       &
                               - 2.0*ec_data(ir  ,1   )       &
                                   + ec_data(ir-1,1   ) ) 
    ec_z2(ir,1 ) = 0.5*dzi**2*(  ec_data(ir  ,2   )         &
                               - ec_data(ir  ,1   ) ) 
  enddo

! Top wall (iz=nz_c): Reflection in z.
  do ir = 1, nr_c-1
    ec_r1(ir,nz_c ) = 0.5*dri*(  ec_data(ir+1,nz_c )          &
                               - ec_data(ir-1,nz_c ) ) 
    ec_z1(ir,nz_c ) = 0.5*dzi*(  ec_data(ir  ,nz_c  )      &
                               - ec_data(ir  ,nz_c-1) )
    ec_r2(ir,nz_c ) = 0.5*dri**2*(      ec_data(ir+1,nz_c)    &
                                  - 2.0*ec_data(ir  ,nz_c)    &
                                      + ec_data(ir-1,nz_c) ) 
    ec_z2(ir,nz_c ) = 0.5*dzi**2*(- ec_data(ir  ,nz_c  )   &
                                  + ec_data(ir  ,nz_c-1) ) 
  enddo

! Corners.
  ec_r1(0,1) =     dri*(  ec_data(1  , 1) )
  ec_z1(0,1) = 0.5*dzi*(  ec_data(1  , 2)    &
                        - ec_data(1  , 1) )
  ec_r2(0,1) = 0.0
  ec_z2(0,1) = 0.5*dzi**2*(  ec_data(1  ,2  )   &
                           - ec_data(1  ,1  ) ) 


  ec_r1(0,nz_c) =     dri*(  ec_data(1  ,nz_c) )
  ec_z1(0,nz_c) = 0.5*dzi*(  ec_data(1  ,nz_c  )    &
                           - ec_data(1  ,nz_c-1) )
  ec_r2(0,nz_c) = 0.0
  ec_z2(0,nz_c) = 0.5*dzi**2*(- ec_data(1  ,nz_c  )   &
                              + ec_data(1  ,nz_c-1) ) 


  ec_r1(nr_c,1) = dri*(- ec_data(nr_c-1,  1) )
  ec_z1(nr_c,1) = 0.5*dzi*(  ec_data(nr_c, 2)    &
                           - ec_data(nr_c, 1) )
  ec_r2(nr_c,1) = 0.0
  ec_z2(nr_c,1) = 0.5*dzi**2*(  ec_data(nr_c,2  )   &
                              - ec_data(nr_c,1  ) ) 



  ec_r1(nr_c,nz_c) = dri*(- ec_data(nr_c-1,nz_c) )
  ec_z1(nr_c,nz_c) = 0.5*dzi*(  ec_data(nr_c,nz_c  )    &
                              - ec_data(nr_c,nz_c-1) )
  ec_r2(nr_c,nz_c) = 0.0
  ec_z2(nr_c,nz_c) = 0.5*dzi**2*(- ec_data(nr_c,nz_c  )   &
                                 + ec_data(nr_c,nz_c-1) ) 

  return

end subroutine sim_rzInterpolateEC
