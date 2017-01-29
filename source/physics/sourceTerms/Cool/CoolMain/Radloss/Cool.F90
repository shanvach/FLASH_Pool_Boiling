!!****if* source/physics/sourceTerms/Cool/CoolMain/Radloss/Cool
!!
!! NAME
!!
!!  Cool(integer(IN)::blockCount
!!       integer(IN)::blockList(blockCount),
!!         real(IN) :: dt,
!!         real(IN) :: time)
!!
!!
!!
!! DESCRIPTION
!!  Apply the radloss cooling opperator on the list of blocks provided as input
!!  The radiative losses rate is used to update the
!!  internal energy in the zone. The radiative losses from an 
!!  optically thin plasma are from Sarazin (1986) Rev.Mod.Phys.
!!  and from Raymond (1976).
!!  (see also Peres et al. 1982, ApJ 252, 791)
!!
!!  After we call radloss, call the eos to update the pressure 
!!  and temperature based on the radiative losses.
!!
!!
!!
!! ARGUMENTS
!!  blockList(:) : The list of blocks on which to apply the cooling operator
!!  blockCount : The number of blocks in the list
!!  dt : the current timestep
!!  time : the current time
!!
!!***


subroutine Cool(blockCount,blockList,dt, time)

!=======================================================================
  use Cool_data, ONLY : cl_Xin,cl_Abar,cl_tradmin,cl_tradmax,&
       cl_dradmin,cl_dradmax,useCool, &
       cl_h1SpecIndex, cl_elecSpecIndex
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_getBlkIndexLimits, &
    Grid_releaseBlkPtr
  use Eos_interface, ONLY : Eos_wrapped
  use Driver_interface, ONLY : Driver_abortFlash
#include "constants.h"
#include "Eos.h"
#include "Flash.h"

  implicit none

  integer, intent(IN) :: blockCount
  integer,dimension(blockCount), intent(IN) :: blockList
  real, intent(IN) :: dt, time
  integer :: blockID
  real, pointer, dimension(:,:,:,:) :: solnData
  integer, dimension(2,MDIM)::blkLimitsGC,blkLimits
  integer :: i, j, k, n
  
  real :: radia,sdot
  real :: hmass
  
  real :: tmp, rho, den_H, den_e, ei, ek
  
  logical :: rad_zone
  integer :: status
  
  if (.not.useCool) return


  do blockID = 1, blockCount


!==============================================================================
     call Grid_getBlkPtr(blockList(blockID),solnData)
     call Grid_getBlkIndexLimits(blockList(blockID),blkLimits,blkLimitsGC)

     call PhysicalConstants_get("proton mass", hmass)

     rad_zone = .FALSE.

! sweep over all the zones

     do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
        do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
           do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              !
              tmp = solnData(TEMP_VAR,i,j,k)
              rho = solnData(DENS_VAR,i,j,k)

! load the mass fractions
              cl_Xin = solnData(SPECIES_BEGIN:SPECIES_BEGIN+NSPECIES-1,i,j,k)
           
           ! derive the Hydrogen and the electron number density
              den_H = rho*cl_Xin(cl_h1SpecIndex-SPECIES_BEGIN+1) &
                       / (cl_Abar(cl_h1SpecIndex-SPECIES_BEGIN+1)*hmass)
              !! DEV : Flash's general Species list does not
              !!       have an entry to electron, the simulation
              !!       directory will have to provide a custom
              !!       implementation of Simulation_initComposition
              if (cl_elecSpecIndex .GE. SPECIES_BEGIN) then
                 den_e = rho*cl_Xin(cl_elecSpecIndex-SPECIES_BEGIN+1) &
                          / (cl_Abar(cl_elecSpecIndex-SPECIES_BEGIN+1)*hmass)
              else
                 call Driver_abortFlash('Cool: electron species not defined!')
              end if

              sdot = 0.0e0
              
              ! if the temperature is between the limits and
              ! the density is between the limits
              
              if ( (tmp >= cl_tradmin .AND. tmp <= cl_tradmax) .AND.&
                   (den_H >= cl_dradmin .AND. den_H <= cl_dradmax) ) then
                 
                 rad_zone = .TRUE.
                 
                 ! radiative losses from an optically thin plasma
                 call radloss(tmp,radia)
                 sdot = -(den_H*den_e*radia)/rho
                 
                 ! change in internal energy due to radiative losses
                 ek = 0.5e0*(solnData(VELX_VAR,i,j,k)**2 + &
                      solnData(VELY_VAR,i,j,k)**2 + &
                      solnData(VELZ_VAR,i,j,k)**2)
                 
                 ei = solnData(EINT_VAR,i,j,k)
                 ei = ei + dt*sdot
                 
                 ! update the global thermodynamic quantities due to the radiative losses
                 solnData(ENER_VAR,i,j,k) = ei + ek
                 solnData(EINT_VAR,i,j,k) = ei
                 
                 ! store the radiative losses rate -- we store it along with the
                 ! solution variables since it is useful to be able to refine on enuc.
                 ! This is the only storage that can currently be refined on.
                 solnData(ENUC_VAR,i,j,k) = sdot
                 
              endif
              
           enddo
        enddo
     enddo
     
     ! if we called the RADLOSS on any zones in this block, then crank the
     ! eos out on the entire block
     if (rad_zone) then 
        call Eos_wrapped(MODE_DENS_EI,blkLimits,blockList(blockID))
     end if
     call Grid_releaseBlkPtr(blockList(blockID),solnData)
  end do
  return
  
end subroutine Cool
