!!****if* source/Simulation/SimulationMain/magnetoHD/CannonballGalaxy/Cool
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
  use Cool_data, ONLY: cl_fctn_file, cl_Nmax,cl_logT, &
       cl_logLambda, cl_Zbar, &
       cl_Nt, cl_na, cl_Xin, cl_Abar, useCool,&
       cl_rho, cl_gamma, cl_smalle, cl_Boltzmann, cl_AMU, cl_Z

  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_getBlkIndexLimits, &
    Grid_releaseBlkPtr, Grid_getCellCoords
  use Eos_interface, ONLY : Eos_wrapped, Eos
  use Multispecies_interface, ONLY : Multispecies_getSumInv, &
       Multispecies_getSumFrac

  use Simulation_data, ONLY: sim_xCtr, &
       sim_yCtr, sim_zCtr, sim_useMeKaLCooling

#include "Eos.h"
#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"

  implicit none
  
  integer,intent(IN)  :: blockCount
  integer,dimension(blockCount), intent(IN) :: blockList
  real,intent(IN) :: dt, time

  real, dimension(EOS_NUM) :: eosData
  logical,dimension(EOS_VARS+1:EOS_NUM) :: eosMask

  integer :: vecLen

  real,pointer,dimension(:,:,:,:) :: solnData

  integer, dimension(LOW:HIGH,MDIM) :: eosRange

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: blockID
  real            :: dtmin, dtguess, t, abund
  real            :: eid, eidold, ei, ek, dedt, r
  integer         :: i, j, k, nok, nbad,lb, size(3)
  
  real, dimension(:), allocatable :: xc, yc, zc
  
  logical, save :: gcell = .true.

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
     size(1) = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
     size(2) = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
     size(3) = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
     allocate(xc(size(1)), yc(size(2)), zc(size(3)))
     call Grid_getCellCoords(KAXIS, blockId, CENTER, gcell, zc, size(3))
     call Grid_getCellCoords(JAXIS, blockId, CENTER, gcell, yc, size(2))
     call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xc, size(1))

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
              call Multispecies_getSumInv(A,cl_Abar,cl_Xin(:))
              cl_Abar = 1.e0 / cl_Abar
              call Multispecies_getSumFrac(Z,cl_Zbar,cl_Xin(:))
              cl_Zbar = cl_Abar*cl_Zbar
#else
              cl_Xin    = 1.0
              abund     = cl_Z
#endif
              cl_gamma  = solnData(GAME_VAR,i,j,k) 

              if (sim_useMeKaLCooling) then

! Start and end times (to avoid roundoff)
 
                 t1     = time
                 t2     = time + dt
                 dt_alt = t2 - t1
                 t2     = t1 + dt_alt
 
! Update the internal energy density using Adams method.

                 call cool_integrate(eid, abund, time, dt)

! If we want to store the cooling in a variable (like for a unit test or something)
! then do it here
#ifdef COOL_VAR
                 call cool_deriv(e, abund, dedt)
                 
                 solnData(COOL_VAR,i,j,k) = dedt ! in units of erg s^-1 cm^-3
#endif              

! update the global thermodynamic quantities due to the radiative losses
                 ei = max(eid/cl_rho, cl_smalle)
                 solnData(ENER_VAR,i,j,k) = ei + ek
                 solnData(EINT_VAR,i,j,k) = ei

              endif

           enddo
        enddo
     enddo

     deallocate(xc)
     deallocate(yc)
     deallocate(zc)
     
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


subroutine cool_deriv (e, abund, dedt)

!===============================================================================

  use Cool_data, ONLY : cl_rho, cl_smalle, cl_logT, cl_abund, &
       cl_gamma, cl_Boltzmann, cl_Abar, cl_logLambda,&
       cl_Nt, cl_Na, cl_da, cl_dlogT, cl_Zbar, cl_AMU, cl_X

  implicit none
  
  real, intent(IN)       :: e, abund
  real, intent(INOUT)    :: dedt
  
  real :: logT, logLambda
  
  integer            :: j, k
  real               :: ne, nH, err, ni, dia, dit,t

  !------------------------------------------------------------------------------
  
  ! Don't cool below the "small energy" threshold.
  
  if (e < cl_rho*cl_smalle) then
     dedt = 0.
     return
  endif
  
  ! OK, find our spot in the table.
  
  t = cl_Abar*(cl_gamma-1.)*cl_AMU*(e/cl_rho)/cl_boltzmann

  logT = log10(t)

  call ut_hunt (cl_logT, cl_Nt, logT, j)
  call ut_hunt (cl_abund, cl_Na, abund, k)

  ! Compute the cooling rate.
  
  ne = 0.5*(1.+cl_X)*cl_rho/cl_AMU
  nH = cl_X*cl_rho/cl_AMU
  ni = (0.25*(1.-cl_X)+cl_X)*cl_rho/cl_AMU

  ! Sanity check, but shouldn't be able to go under 
  ! zero metallicity anyway

  if (k == 0) k = 1

  if (j == 0) then
     
     ! turn off cooling at low temperatures to avoid runaway
     
     dedt = 0.
     
  elseif (j == cl_Nt) then
     
     ! free-free thermal emission at high temperatures
     ! (Rybicki & Lightman 1979 eq. 5.15b)
     
     dedt = -1.4E-27 * sqrt(t) * ne * ni * cl_Zbar * cl_Zbar
     
  else
     
     ! linear interpolation (in log(T) and a) from table
     
     dit = (logT-cl_logT(j))/cl_dlogT
     dia = (abund-cl_abund(k))/cl_da

     logLambda = cl_logLambda(j,k)*(1.-dit)*(1.-dia) + &
          cl_logLambda(j+1,k)*(1.-dia)*dit + &
          cl_logLambda(j,k+1)*dia*(1.-dit) + &
          cl_logLambda(j+1,k+1)*dit*dia

     dedt = - (ne*nH) * (10.**logLambda)
     
  endif
  
  !===============================================================================
  
  return
end subroutine cool_deriv



