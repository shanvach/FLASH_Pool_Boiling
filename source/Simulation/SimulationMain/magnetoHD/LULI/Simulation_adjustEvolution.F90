!!****if* source/Simulation/SimulationMain/magnetoHD/LULI/Simulation_adjustEvolution
!!
!! NAME
!!  Simulation_adjustEvolution
!!
!!
!! SYNOPSIS
!!  Simulation_adjustEvolution( integer(IN) :: blkcnt,
!!                              integer(IN) :: blklst(blkcnt),
!!                              integer(IN) :: nstep,
!!                              real(IN) :: dt,
!!                              real(IN) :: stime )
!!
!! DESCRIPTION
!!  This routine is called every cycle. It can be used to adjust
!!  the simulation while it is running.
!!  
!! ARGUMENTS
!!  blkcnt - number of blocks
!!  blklist - block list
!!  nstep - current cycle number
!!  dt - current time step length
!!  stime - current simulation time
!!
!!***
subroutine Simulation_adjustEvolution(blkcnt, blklst, nstep, dt, stime)
  use Simulation_data, ONLY : sim_targetMass, sim_meshComm, sim_pulseLength, &
                              sim_laserEnergy, sim_computeBiermann, &
                              sim_speedlt, sim_qele, sim_avo, sim_driverType
  use Grid_interface,  ONLY : Grid_getSingleCellVol, &
                              Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_fillGuardCells, &
                              Grid_releaseBlkPtr, Grid_getSingleCellVol
  use Eos_interface,   ONLY : Eos_wrapped, Eos_getAbarZbar
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer, intent(in) :: blkcnt
  integer, intent(in) :: blklst(blkcnt)
  integer, intent(in) :: nstep
  real, intent(in) :: dt
  real, intent(in) :: stime
  integer :: blkLimitsGC(LOW:HIGH,MDIM)
  integer :: blkLimits(LOW:HIGH,MDIM)
  integer :: i,j,k, lb
  real, pointer   :: blkPtr(:,:,:,:)
  real :: vol
  real :: targ_mass, mass_frac
  real :: rho
  real :: eele
  real :: elas
  real :: mfrac
  real :: abar, zbar
  real :: gradpe(MDIM), gradne(MDIM), cross(MDIM), bier(MDIM)
  real :: delta(MDIM) ! The cell width

  real :: block_mass, target_mass

  integer :: ierr
  integer :: is, ie, js, je, ks, ke
  real, allocatable :: nele(:,:,:)

  ! *************************************
  ! *                                   *
  ! *     COMPUTE TOTAL TARGET MASS     *
  ! *                                   *
  ! *************************************
  sim_targetMass = 0.0
  target_mass = 0.0

  do lb = 1, blkcnt
     call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blklst(lb), blkPtr)

     block_mass = 0.0
     
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              call Grid_getSingleCellVol(blklst(lb), EXTERIOR, (/i,j,k/), vol)
              mfrac = blkPtr(HEAT_MSCALAR,i,j,k)
              rho = blkPtr(DENS_VAR,i,j,k)              
              block_mass = block_mass + vol * rho * mfrac
           enddo
        end do
     end do
     call Grid_releaseBlkPtr(blklst(lb), blkPtr)

     target_mass = target_mass + block_mass
  end do  

  ! Sum the target mass across all processes:  
  call mpi_allreduce(target_mass, sim_targetMass, 1, FLASH_REAL, & 
       MPI_SUM, sim_meshComm, ierr )     

  ! *******************************************
  ! *                                         *
  ! *     COMPUTE BIERMANN BATTERY SOURCE     *
  ! *                                         *
  ! *******************************************
  if(sim_computeBiermann) then

     call Grid_fillGuardCells(CENTER,ALLDIR)

     do lb = 1, blkcnt
        call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
        call Grid_getBlkPtr(blklst(lb), blkPtr)
        call Grid_getDeltas(blklst(lb), delta)

        is = blkLimitsGC(LOW,  IAXIS)
        js = blkLimitsGC(LOW,  JAXIS)
        ks = blkLimitsGC(LOW,  KAXIS)
        ie = blkLimitsGC(HIGH, IAXIS)
        je = blkLimitsGC(HIGH, JAXIS)
        ke = blkLimitsGC(HIGH, KAXIS)
        
        allocate(nele(is:ie, js:je, ks:ke))
        
        ! Compute nele in each cell for this block:
        do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
           do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
              do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
                 call Eos_getAbarZbar(blkPtr(:,i,j,k), abar=abar, zbar=zbar)
                 nele(i,j,k) = zbar/abar * sim_avo * blkPtr(DENS_VAR,i,j,k)
              enddo
           end do
        end do

        ! Compute the Biermann Battery source:
        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                 gradpe    = 0.0
                 gradne    = 0.0
                 cross     = 0.0

                 ! Compute gradient of electron pressure and gradient of
                 ! electron number density:
                 gradpe(IAXIS) = (blkPtr(PELE_VAR,i+1,j,k) - blkPtr(PELE_VAR,i-1,j,k))/(2*delta(IAXIS))
                 gradne(IAXIS) = (nele(i+1,j,k) - nele(i-1,j,k))/(2*delta(IAXIS))

#if NDIM == 2 || NDIM == 3
                 gradpe(JAXIS) = (blkPtr(PELE_VAR,i,j+1,k) - blkPtr(PELE_VAR,i,j-1,k))/(2*delta(JAXIS))
                 gradne(JAXIS) = (nele(i,j+1,k) - nele(i,j-1,k))/(2*delta(JAXIS))
#endif

#if NDIM == 3
                 gradpe(KAXIS) = (blkPtr(PELE_VAR,i,j,k+1) - blkPtr(PELE_VAR,i,j,k-1))/(2*delta(KAXIS))
                 gradne(KAXIS) = (nele(i,j,k+1) - nele(i,j,k-1))/(2*delta(KAXIS))
#endif

                 ! Compute grad(PELE) cross grad(NELE):
#if NDIM == 2 || NDIM == 3
                 cross(KAXIS) = gradpe(IAXIS)*gradne(JAXIS) - gradpe(JAXIS)*gradne(IAXIS)
#endif

#if NDIM == 3
                 cross(IAXIS) = gradpe(JAXIS)*gradne(KAXIS) - gradpe(KAXIS)*gradne(JAXIS)
                 cross(JAXIS) = gradpe(KAXIS)*gradne(IAXIS) - gradpe(IAXIS)*gradne(KAXIS)
#endif

                 ! Compute the Biermann Battery term:
                 ! dB/dt = c/nele**2/e * (grad(pele) cross grad(nele))
                 bier = sim_speedlt * cross / (nele(i,j,k)**2 * sim_qele)
                 blkPtr(BIEX_VAR,i,j,k) = bier(IAXIS)
                 blkPtr(BIEY_VAR,i,j,k) = bier(JAXIS)
                 blkPtr(BIEZ_VAR,i,j,k) = bier(KAXIS)
                 
                 blkPtr(GNEX_VAR,i,j,k) = gradne(IAXIS)
                 blkPtr(GNEY_VAR,i,j,k) = gradne(JAXIS)
                 blkPtr(GNEZ_VAR,i,j,k) = gradne(KAXIS)

                 blkPtr(GPEX_VAR,i,j,k) = gradpe(IAXIS)
                 blkPtr(GPEY_VAR,i,j,k) = gradpe(JAXIS)
                 blkPtr(GPEZ_VAR,i,j,k) = gradpe(KAXIS)

                 blkPtr(NELE_VAR,i,j,k) = nele(i,j,k)


              enddo
           end do
        end do

        deallocate(nele)

        call Grid_releaseBlkPtr(blklst(lb), blkPtr)
     end do
  end if

  ! Deposit energy:

  if(stime < sim_pulseLength .and. &
       (sim_driverType == "uniform" .or. sim_driverType == "unispec")) then

     ! Compute amount of energy to dump in on this time
     ! step:
     elas = sim_laserEnergy * dt/sim_pulseLength

     if (sim_driverType == "uniform") then
        elas = elas/sim_targetMass
     end if

     do lb = 1, blkcnt
        call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
        call Grid_getBlkPtr(blklst(lb), blkPtr)
        
        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
                 
                 
                 eele  = blkPtr(EELE_VAR, i,j,k) 
                 mfrac = blkPtr(HEAT_MSCALAR,i,j,k)
                 
                 eele = eele + mfrac*elas
                 blkPtr(EELE_VAR,i,j,k) = eele
              enddo
           end do
        end do
        call Grid_releaseBlkPtr(blklst(lb), blkPtr)
        
        call Eos_wrapped(MODE_DENS_EI_GATHER,blkLimits,blklst(lb))
     end do
  end if

end subroutine Simulation_adjustEvolution
