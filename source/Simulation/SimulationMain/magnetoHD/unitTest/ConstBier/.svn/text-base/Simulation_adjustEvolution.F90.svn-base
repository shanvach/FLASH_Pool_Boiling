!!****if* source/Simulation/SimulationMain/magnetoHD/unitTest/ConstBier/Simulation_adjustEvolution
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
  use Simulation_data, ONLY : sim_nele1, sim_nele2, sim_pele1, sim_pele2, sim_xmin, &
                              sim_xmax, sim_ymin, sim_ymax, sim_avogadro, sim_qele, &
                              sim_speedlt
  use Grid_interface,  ONLY : Grid_getSingleCellVol, &
                              Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_fillGuardCells, &
                              Grid_releaseBlkPtr
  use Eos_interface,   ONLY : Eos_wrapped, Eos_getAbarZbar
  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: blkcnt
  integer, intent(in) :: blklst(blkcnt)
  integer, intent(in) :: nstep
  real, intent(in) :: dt
  real, intent(in) :: stime
  integer :: blkLimitsGC(LOW:HIGH,MDIM)
  integer :: blkLimits(LOW:HIGH,MDIM)
  integer :: i,j,k, lb
  real, pointer   :: blkPtr(:,:,:,:)
  integer :: point(MDIM)
  real :: vol
  real :: targ_mass, mass_frac
  real :: rho
  real :: eele
  real :: elas
  real :: mfrac
  real :: abar, zbar
  real :: delta(MDIM) ! The cell width

  real :: gradne
  real :: gradpe
  real :: nele

  ! call Grid_fillGuardCells(CENTER,ALLDIR)

  gradne = (sim_nele2-sim_nele1)/(sim_ymax-sim_ymin)
  gradpe = (sim_pele2-sim_pele1)/(sim_xmax-sim_xmin)

  do lb = 1, blkcnt
     call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blklst(lb), blkPtr)
     
     ! Compute the Biermann Battery source:
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
              call Eos_getAbarZbar(blkPtr(:,i,j,k), abar=abar, zbar=zbar)
              nele = zbar*sim_avogadro*blkPtr(DENS_VAR,i,j,k)/abar
              blkPtr(BIER_VAR,i,j,k) = sim_speedlt * (stime + dt) * gradne * gradpe / (nele**2 * sim_qele)

           enddo
        end do
     end do
     
     call Grid_releaseBlkPtr(blklst(lb), blkPtr)
  end do

end subroutine Simulation_adjustEvolution
