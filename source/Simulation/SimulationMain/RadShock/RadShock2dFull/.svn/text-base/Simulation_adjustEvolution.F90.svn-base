!!****if* source/Simulation/SimulationMain/RadShock/RadShock2dFull/Simulation_adjustEvolution
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
  use rt_data, ONLY: rt_mgdDomainBC, rt_useMGD
  use Logfile_interface, ONLY: Logfile_stampMessage
  implicit none

#include "constants.h"

  integer, intent(in) :: blkcnt
  integer, intent(in) :: blklst(blkcnt)
  integer, intent(in) :: nstep
  real, intent(in) :: dt
  real, intent(in) :: stime

  ! Turn off the X-ray drive after 1 ns:
  if(rt_useMGD .eqv. .true.) then
     if(stime > 1.0e-09 .and. any(rt_mgdDomainBC(:,3) /= VACUUM) ) then
        call Logfile_stampMessage("Deactivating X-Ray Drive")
        rt_mgdDomainBC(:,3) = VACUUM
     end if
  end if


end subroutine Simulation_adjustEvolution
