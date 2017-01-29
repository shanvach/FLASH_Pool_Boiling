!!****if* source/Simulation/SimulationMain/magnetoHD/AlfvenWave/hy_uhd_energyFix
!!
!! NAME
!!
!!  hy_uhd_energyFix
!!
!! SYNOPSIS
!!
!!  hy_uhd_energyFix( integer (IN) :: blockID,
!!                    integer (IN) :: blkLimits(2,MDIM),
!!                    real(IN)     :: dt,
!!                    real(IN)     :: del(MDIM),
!!                    integer(IN)  :: eosMode)
!!
!! DESCRIPTION
!!
!!  This routine corrects energy in two different ways:
!!  The first choice is to fix the energy due to the differences of
!!  the magnetic pressures using the cell-centered magnetic 
!!  fields and the divergence-free cell face-centered magnetic fields.
!!  This correction is optional, but may be useful for low beta
!!  plasma flows. To enable this first correction during the simulation
!!  two runtime parameters "hy_killdivb" and "hy_energyFixSwitch" 
!!  should be both turned on in flash.par file.
!!  The second correction is to use the internal energy evolution
!!  to avoid any negativity states of pressure in the Eos routines.
!!
!!  Note that this routine does more than just correcting energy, and
!!  also computes several quantities such as divergence of magnetic
!!  fields, total pressure, current density, and electric fields, etc.
!!  This routine also takes care of abundances.
!!
!! ARGUMENTS
!!
!!  blockID   - a local block ID
!!  blkLimits - an array that holds the lower and upper indices of the section
!!              of block without the guard cells
!!  dt        - time step
!!  del       - grid deltas in each direction
!!  eosMode   - a mode used in a call to Eos
!!
!!***

!!REORDER(4):E,U,B[xyz]

Subroutine hy_uhd_energyFix(blockID,blkLimits,dt,del,eosMode)

  use Hydro_data,     ONLY : hy_killdivb, hy_eswitch, &
                             hy_energyFixSwitch, hy_irenorm

  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr,&
                             Grid_renormAbundance, Grid_limitAbundance                             

  implicit none

#include "Flash.h"
#include "constants.h"
#include "UHD.h"

  !! ---- Argument List ----------------------------------
  integer, intent(IN) :: blockID
  integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits
  real, intent(IN) :: dt
  real, dimension(MDIM), intent(IN) :: del
  integer, intent(IN) :: eosMode
  !! -----------------------------------------------------

  integer :: i,j,k
  real    :: ekin,eint,emag, mhdEnergyCorrection
  real, pointer, dimension(:,:,:,:) :: U
  real :: cos_ang, sin_ang
#if NFACE_VARS > 0
#if NDIM > 1
  real, pointer, dimension(:,:,:,:) :: E,Bx,By,Bz
#endif
#endif

!!$  real :: ye,nele_left,nele_rght,nele,pele_left,pele_rght,pele
!!$  real :: gradNelex,gradPelex,gradNeley,gradPeley


  call Grid_getBlkPtr(blockID,U,CENTER)
#if NFACE_VARS > 0
#if NDIM > 1
  call Grid_getBlkPtr(blockID,E,SCRATCH)
  call Grid_getBlkPtr(blockID,Bx,FACEX)
  call Grid_getBlkPtr(blockID,By,FACEY)
  if (NDIM == 3) call Grid_getBlkPtr(blockID,Bz,FACEZ)
#endif
#endif

  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

#ifdef BDRY_VAR
           if (U(BDRY_VAR,i,j,k) .ne. 1.) then
#endif
           ekin = .5*dot_product(U(VELX_VAR:VELZ_VAR,i,j,k),U(VELX_VAR:VELZ_VAR,i,j,k))&
                    *U(DENS_VAR,i,j,k)
#ifdef BDRY_VAR
           if (U(BDRY_VAR,i,j,k) == 1.) ekin = 0.
#endif
           !! emag is a magnetic pressure calculated using divergence-free
           !! face-centered B fields on the staggered grid and take 
           !! an arithmetic average to get cell-centered B fields 
           !! (see hy_uhd_staggeredDivb.F90)

           !! MAGP_VAR is a magnetic pressure calculated using a standard Godunov
           !! update based on a cell-centered B fields, which are non-divergence-free.

           emag = .5*dot_product(U(MAGX_VAR:MAGZ_VAR,i,j,k),U(MAGX_VAR:MAGZ_VAR,i,j,k))

           if (hy_killdivb .and. hy_energyFixSwitch) then
              !! U(MAGP_VAR,i,j,k) is a magnetic pressure from the Godunov cell-centered update
              mhdEnergyCorrection = emag - U(MAGP_VAR,i,j,k)
              U(ENER_VAR,i,j,k) = U(ENER_VAR,i,j,k)+mhdEnergyCorrection
           endif

           U(MAGP_VAR,i,j,k) = emag
           eint = U(ENER_VAR,i,j,k)-ekin-emag


           !! (1) if epsilon --> 0, then we use the total conservative energy to get eint;
           !! (2) if epsilon >> 0, then the eint is obtained by solving a non-conservative 
           !!     auxilary internal energy equation.
           if ((eint >= hy_eswitch*(ekin+emag) ) .and. (eint >= hy_eswitch*U(EINT_VAR,i,j,k))) then
              U(ENER_VAR,i,j,k) = (ekin + eint)/U(DENS_VAR,i,j,k)
              U(EINT_VAR,i,j,k) = eint/U(DENS_VAR,i,j,k)
           else
              !! Store specific gas energy, not plasma energy
              U(ENER_VAR,i,j,k) = (ekin + U(EINT_VAR,i,j,k))/U(DENS_VAR,i,j,k)
              !! Specific internal energy
              U(EINT_VAR,i,j,k) = U(EINT_VAR,i,j,k)/U(DENS_VAR,i,j,k)
           endif

#ifdef DIVV_VAR
#if NDIM == 1
           U(DIVV_VAR,i,j,k) = (U(VELX_VAR,i+1,j,k)-U(VELX_VAR,i-1,j,k))/del(DIR_X)
#elif NDIM == 2
           U(DIVV_VAR,i,j,k) = (U(VELX_VAR,i+1,j,k)-U(VELX_VAR,i-1,j,k))/del(DIR_X)&
                              +(U(VELY_VAR,i,j+1,k)-U(VELY_VAR,i,j-1,k))/del(DIR_Y)
#elif NDIM == 3
           U(DIVV_VAR,i,j,k) = (U(VELX_VAR,i+1,j,k)-U(VELX_VAR,i-1,j,k))/del(DIR_X)&
                              +(U(VELY_VAR,i,j+1,k)-U(VELY_VAR,i,j-1,k))/del(DIR_Y)&
                              +(U(VELZ_VAR,i,j,k+1)-U(VELZ_VAR,i,j,k-1))/del(DIR_Z)
#endif
#endif


#if NFACE_VARS > 0
#if NDIM == 1
           U(DIVB_VAR,i,j,k) =  0.
#elif NDIM == 2
           U(DIVB_VAR,i,j,k) = &
                (Bx(MAG_FACE_VAR,i+1,j,  k) - Bx(MAG_FACE_VAR,i,j,k))/del(DIR_X) &
              + (By(MAG_FACE_VAR,i,  j+1,k) - By(MAG_FACE_VAR,i,j,k))/del(DIR_Y)
#elif NDIM == 3
           U(DIVB_VAR,i,j,k) = &
                (Bx(MAG_FACE_VAR,i+1,j,  k  ) - Bx(MAG_FACE_VAR,i,j,k))/del(DIR_X) &
              + (By(MAG_FACE_VAR,i,  j+1,k  ) - By(MAG_FACE_VAR,i,j,k))/del(DIR_Y) &
              + (Bz(MAG_FACE_VAR,i,  j,  k+1) - Bz(MAG_FACE_VAR,i,j,k))/del(DIR_Z)
#endif
#else 
           !! pure hydro mode
           U(DIVB_VAR,i,j,k) = 0.
#endif

#ifdef TOTP_VAR
           !! Total plasma pressure and plasma beta
           U(TOTP_VAR,i,j,k) = U(PRES_VAR,i,j,k)+U(MAGP_VAR,i,j,k)
#endif

#ifdef BETA_VAR
           U(BETA_VAR,i,j,k) = U(PRES_VAR,i,j,k)/U(MAGP_VAR,i,j,k)
#endif

#if NFACE_VARS > 0
#if NDIM > 1
#ifdef VECZ_VAR
           !! advance Az to next time step in 2D
           !! dA/dt = -E = VxB - magVisc*J
           U(VECZ_VAR,i,j,k) = U(VECZ_VAR,i,j,k)&
                -0.25*dt*( E(EZ_SCRATCH_GRID_VAR, i,  j,  k) &
                          +E(EZ_SCRATCH_GRID_VAR, i+1,j,  k) &
                          +E(EZ_SCRATCH_GRID_VAR, i,  j+1,k) &
                          +E(EZ_SCRATCH_GRID_VAR, i+1,j+1,k))
#endif
#ifdef ELEX_VAR
           U(ELEX_VAR,i,j,k) =( E(EX_SCRATCH_GRID_VAR, i,j,  k  ) &
                               +E(EX_SCRATCH_GRID_VAR, i,j,  k+1) &
                               +E(EX_SCRATCH_GRID_VAR, i,j+1,k  ) &
                               +E(EX_SCRATCH_GRID_VAR, i,j+1,k+1))*0.25
#endif
#ifdef ELEY_VAR
           U(ELEY_VAR,i,j,k) =( E(EY_SCRATCH_GRID_VAR, i,  j,k  ) &
                               +E(EY_SCRATCH_GRID_VAR, i+1,j,k  ) &
                               +E(EY_SCRATCH_GRID_VAR, i,  j,k+1) &
                               +E(EY_SCRATCH_GRID_VAR, i+1,j,k+1))*0.25
#endif
#ifdef ELEZ_VAR
           U(ELEZ_VAR,i,j,k) =( E(EZ_SCRATCH_GRID_VAR, i,  j,  k) &
                               +E(EZ_SCRATCH_GRID_VAR, i+1,j,  k) &
                               +E(EZ_SCRATCH_GRID_VAR, i,  j+1,k) &
                               +E(EZ_SCRATCH_GRID_VAR, i+1,j+1,k))*0.25
#endif
#endif
#endif


!!$#ifdef XNOR_VAR
!!$           solnData(XNOR_VAR,i,j,k) = x1
!!$#endif
#ifdef BTRA_VAR
           sin_ang=2./sqrt(5.)
           cos_ang=sin_ang*0.5
           U(BTRA_VAR,i,j,k) =-U(MAGX_VAR,i,j,k)*sin_ang+U(MAGY_VAR,i,j,k)*cos_ang
#endif


#ifdef BDRY_VAR
        endif
#endif



       enddo
     enddo
  enddo


  if (eosMode==MODE_DENS_PRES) then
     U(PRES_VAR,:,:,:) = U(EINT_VAR,:,:,:)*U(DENS_VAR,:,:,:)*(U(GAME_VAR,:,:,:)-1.)
  endif

  ! Renormalize or limit abundances
  if (hy_irenorm == 1) then
     call Grid_renormAbundance(blockID,blkLimits,U)
  else
     call Grid_limitAbundance(blkLimits,U)
  endif
  

  call Grid_releaseBlkPtr(blockID,U,CENTER)
#if NFACE_VARS > 0
#if NDIM > 1
  call Grid_releaseBlkPtr(blockID,E,SCRATCH)
  call Grid_releaseBlkPtr(blockID,Bx,FACEX)
  call Grid_releaseBlkPtr(blockID,By,FACEY)
  if (NDIM == 3) call Grid_releaseBlkPtr(blockID,Bz,FACEZ)
#endif
#endif

  return
End Subroutine hy_uhd_energyFix
