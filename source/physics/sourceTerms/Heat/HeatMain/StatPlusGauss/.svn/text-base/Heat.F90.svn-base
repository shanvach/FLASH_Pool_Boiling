!!****if* source/physics/sourceTerms/Heat/HeatMain/StatPlusGauss/Heat
!!
!! NAME
!!  
!!  Heat 
!!
!!
!! SYNOPSIS
!! 
!!  call Heat (integer(IN) :: blockCount,
!!             integer(IN) :: blockList(blockCount),
!!             real(IN)    :: dt,
!!             real(IN)    :: time)
!!
!!  
!!  
!! DESCRIPTION
!!
!!  Apply the stat+gauss source term operator to a block of
!!  zones. The energy generation rate is used to update the
!!  internal energy in the zone. The phenomenological heating
!!  rate is described as a 3-D Gauss function.
!!
!!  After we call stat+gauss, call the eos to update the
!!  pressure and temperature based on the phenomenological
!!  heating.
!!  
!!
!! ARGUMENTS
!!
!!  blockCount : number of blocks to operate on
!!  blockList  : list of blocks to operate on
!!  dt         : current timestep
!!  time       : current time
!!
!!***
subroutine Heat (blockCount,blockList,dt,time)
!
!==============================================================================
!!$  use Heat_data, ONLY: ht_x0, ht_y0, ht_z0,&
!!$       ht_stat, ht_q, ht_sig, ht_tstar, ht_t0, ht_tau, ht_tmin, &
!!$       useHeat, ht_meshMe, ht_numProcs

  use Heat_data
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_releaseBlkPtr
  use Eos_interface, ONLY : Eos_wrapped
!
#include "constants.h"
#include "Flash.h"
  implicit none
  
  integer,intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN)::blockList
  real,intent(IN) :: dt,time
  integer :: blockID,lb
  real,pointer,dimension(:,:,:,:):: solnData
  real :: tranheat, argm, sdot
  !
  real :: tmp, rho, ei, ek
  !
  real,allocatable, dimension(:) :: xCenter, yCenter, zCenter
  real :: xx, yy, zz
  integer,dimension(MDIM)  :: dimSize
  logical :: gcell = .true.
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC

  integer :: i,j,k
  !==============================================================================
!
  if(.not.useHeat) return


  do lb = 1, blockCount
     blockID=blockList(lb)
     !
     call Grid_getBlkPtr(blockID,solnData)
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     dimSize(:)=blkLimitsGC(HIGH,:)-blkLimitsGC(LOW,:)+1
     if (NDIM > 2)then
        allocate(zCenter(dimSize(KAXIS)))
        call Grid_getCellCoords(KAXIS,blockID,&
                                 CENTER,gcell,zCenter,dimSize(KAXIS))
     end if
     if (NDIM > 1)then
        allocate(yCenter(dimSize(JAXIS)))
        call Grid_getCellCoords(JAXIS,blockID,&
                                 CENTER,gcell,yCenter,dimSize(JAXIS))
     end if
     
     allocate(xCenter(dimSize(IAXIS)))
     call Grid_getCellCoords(IAXIS,blockID,&
                                 CENTER,gcell,xCenter,dimSize(IAXIS))



     do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
        if (NDIM == 3) zz = zCenter(k)
        !
        do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
           if (NDIM >= 2) yy = yCenter(j)
           
           do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              xx = xCenter(i)
              !
              tmp = solnData(TEMP_VAR,i,j,k)
              rho = solnData(DENS_VAR,i,j,k)
              !
              sdot = 0.0e0
              if (tmp >= ht_tmin) then
                 !
                 ! Phenomenological heating (Gauss function)
                 tranheat = 0.0
                 if (time >= ht_tstar) then 
                    argm = (xx-ht_x0)**2
                    if (NDIM >= 2) argm = argm + (yy-ht_y0)**2
                    if (NDIM == 3) argm = argm + (zz-ht_z0)**2
                    tranheat = ht_q*exp(-0.5*argm/ht_sig**2)
                    if (time >= ht_t0)                                 &
                         &                 tranheat = tranheat*exp((ht_t0-time)/ht_tau)
                 endif
                 !
                 sdot = (ht_stat+tranheat)/rho
              endif
              !
              ! change in internal energy due to phenomenological heating
              ek = 0.5e0*(solnData(VELX_VAR,i,j,k)**2 +                   &
                   &                    solnData(VELY_VAR,i,j,k)**2 +                   &
                   &                    solnData(VELZ_VAR,i,j,k)**2)
              !
              ei = solnData(ENER_VAR,i,j,k) - ek
              ei = ei + dt*sdot
              !
              ! update the global thermodynamic quantities due to the phen. heating
              solnData(ENER_VAR,i,j,k) = ei + ek
#ifdef EINT_VAR
              solnData(EINT_VAR,i,j,k) = ei
#endif
              !
              ! store the phenomenological heating rate -- we store it along with the
              ! solution variables since it is useful to be able to refine on enuc.
              ! This is the only storage that can currently be refined on.
              !               solnData(ENUC_VAR,i,j,k) = sdot
              !
           enddo
        enddo
     enddo
     !
     ! if we called the STAT+GAUSS on any zones in this block, then crank the
     ! eos out on the entire block
     call Eos_wrapped(MODE_DENS_EI,blkLimits,blockID)
     
     call Grid_releaseBlkPtr(blockID,solndata)
     deallocate(xCenter)
     if(NDIM>1)deallocate(yCenter)
     if(NDIM>2)deallocate(zCenter)
  enddo
  
  return
end subroutine Heat


