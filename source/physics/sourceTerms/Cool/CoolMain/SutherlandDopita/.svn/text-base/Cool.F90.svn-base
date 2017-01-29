!!****if* source/physics/sourceTerms/Cool/CoolMain/SutherlandDopita/Cool
!!
!! NAME
!!
!!  Cool
!!
!! SYNOPSIS
!!
!!  Cool(integer,intent(IN)  :: blockcount,
!!       integer,dimension(blockCount), intent(IN)  :: blocklist,
!!       real,intent(IN)  :: dt,
!!       real,intent(IN)  :: time)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   blockcount : 
!!
!!   blocklist : 
!!
!!   dt : 
!!
!!   time : 
!!
!! AUTOGENROBODOC
!!
!!
!!***


subroutine Cool (blockCount, blockList, dt, time)

!==============================================================================
  use Cool_data, ONLY :cl_fctn_file,cl_Nmax,cl_logT, cl_ne, cl_nH, &
       cl_nt, cl_logLnet, cl_logLnorm, &
       cl_logU, cl_logtau, cl_P12, cl_rho24, cl_Ci, cl_mubar, &
       cl_dlogLdlogU, cl_dnedlogU, cl_dntdlogU,cl_N,cl_Xin, cl_Abar, useCool,&
       cl_rho, cl_gamma, cl_smalle, cl_Boltzmann, cl_AMU

  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_getBlkIndexLimits, &
    Grid_releaseBlkPtr
  use Eos_interface, ONLY : Eos_wrapped

#include "constants.h"
#include "Flash.h"

  implicit none
  
  integer,intent(IN)  :: blockCount
  integer,dimension(blockCount), intent(IN) :: blockList
  real,intent(IN) :: dt, time

  real,pointer,dimension(:,:,:,:) :: solnData
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: blockID
  real            :: dtmin, dtguess
  real            :: eid, eidold, ei, ek, dedt
  integer         :: i, j, k, nok, nbad,lb
  
  
  ! Stuff for ODE solve
  real                  :: t1, t2, dt_alt
  
  external cool_deriv

!==============================================================================
  if(.not.useCool) return


  do lb = 1, blockCount
!!DEV: blockList(blockCount) in the next 2 lines looks really wrong - KW
     blockID=blockList(lb)
     call Grid_getBlkPtr(blockID,solnData)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              cl_rho    = solnData(DENS_VAR,i,j,k)
              eid    = cl_rho * solnData(EINT_VAR,i,j,k)
              eidold = eid
              ek     = 0.5*(solnData(VELX_VAR,i,j,k)**2 + &
                   solnData(VELY_VAR,i,j,k)**2 + &
                   solnData(VELZ_VAR,i,j,k)**2)

#if NSPECIES > 0
              cl_Xin    = solnData(SPECIES_BEGIN:SPECIES_BEGIN+NSPECIES-1,i,j,k)
#else
              cl_Xin    = 1.0
#endif
              cl_gamma  = solnData(GAME_VAR,i,j,k) 
              !!should really get from EOS during integ

! Start and end times (to avoid roundoff)
 
              t1     = time
              t2     = time + dt
              dt_alt = t2 - t1
              t2     = t1 + dt_alt
 
! Update the internal energy density using Adams method.

              call cool_integrate(eid, time, dt)

! If we want to store the cooling in a variable (like for a unit test or something)
! then do it here
#ifdef COOL_VAR
              call cool_deriv(eid, dedt)

              solnData(COOL_VAR,i,j,k) = dedt ! in units of erg s^-1 cm^-3
#endif              

! update the global thermodynamic quantities due to the radiative losses
              ei = max(eid/cl_rho, cl_smalle)
              solnData(ENER_VAR,i,j,k) = ei + ek
              solnData(EINT_VAR,i,j,k) = ei

! store the radiative loss rate -- we store it along with the
! solution variables since it is useful to be able to refine on enuc.
#ifdef ENUC_VAR
              solnData(ENUC_VAR,i,j,k) = (eid - eidold) / dt
#endif              
           enddo
        enddo
     enddo
     
     ! Make all thermodynamic quantities consistent with each other.
     
     call Eos_wrapped(MODE_DENS_EI,blkLimits,blockID)
     call Grid_releaseBlkPtr(blockID,solnData)
  end do
  !===============================================================================
  
  return
end subroutine Cool



!******************************************************************************

! Routine:     cool_deriv

! Description: Routine fed to odeint to supply cooling source term.
!              Interpolates from tables supplied through CoolDataModule.


subroutine cool_deriv (e, dedt)

!===============================================================================

  use Cool_data, ONLY : cl_rho, cl_smalle, cl_rho24, cl_logU,&
       cl_ne,cl_nt,cl_gamma, cl_Boltzmann, cl_mubar,cl_logLnorm,&
       cl_N, cl_dnedlogU, cl_dntdlogU, cl_dlogLdlogU
  implicit none
  
  real, intent(INOUT)    :: e, dedt
  
  real    :: tmp, logtmp, logLambda
  
  integer            :: j, k
  integer, parameter :: m = 4
  real               :: logesc, den_e, den_t, err, rho_scaling

  !------------------------------------------------------------------------------
  
  ! Sutherland & Dopita tables are for constant ion number density, varying
  ! temperature.  The average particle mass is itself a function of temperature,
  ! so we can't compute temperature directly.  Instead, we find the factor by
  ! which we need to scale the (mass) density to bring it into the table, then
  ! compute the corresponding scaled energy density.  (We're effectively
  ! ignoring the mass of the electron here.)
  
  ! Don't cool below the "small energy" threshold.
  
  if (e < cl_rho*cl_smalle) then
     dedt = 0.
     return
  endif
  
  ! OK, find our spot in the table.
  
  rho_scaling = cl_rho / (cl_rho24(1)*1.E-24)
  tmp = e / rho_scaling
  logesc = log10(tmp)

  call ut_hunt (cl_logU, cl_N, logesc, j)
  
  ! Compute the cooling rate.
  
  if (j == 0) then
     
     ! turn off cooling at low temperatures to avoid runaway
     
     dedt = 0.
     
  elseif (j == cl_N) then
     
     ! free-free thermal emission at high temperatures
     ! (Rybicki & Lightman 1979 eq. 5.15b)
     
     den_t = cl_rho / (cl_rho24(cl_N)*1.E-24/cl_nt(cl_N))
     den_e = den_t * cl_ne(cl_N)/cl_nt(cl_N)
     tmp   = e*(cl_gamma-1.)*cl_mubar(cl_N)*1.E-24/(cl_rho*cl_Boltzmann)
     dedt  = -1.91E-27 * sqrt(tmp) * den_e * den_t
     
  else
     
     ! linear interpolation (in log(energy)) from Sutherland & Dopita tables
     
     logLambda = cl_logLnorm(j) + cl_dlogLdlogU(j)*(logesc-cl_logU(j))
     den_e     = cl_ne(j)       + cl_dnedlogU(j)  *(logesc-cl_logU(j))
     den_t     = cl_nt(j)       + cl_dntdlogU(j)  *(logesc-cl_logU(j))
     
     den_e = den_e * rho_scaling
     den_t = den_t * rho_scaling
     dedt  = -(den_t*den_e) * (10.**logLambda)
     
  endif
  
  !===============================================================================
  
  return
end subroutine cool_deriv



