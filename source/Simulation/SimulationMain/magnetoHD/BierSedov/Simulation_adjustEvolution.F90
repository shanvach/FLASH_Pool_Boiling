!!****if* source/Simulation/SimulationMain/magnetoHD/BierSedov/Simulation_adjustEvolution
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
!!  blklst - block list
!!  nstep - current cycle number
!!  dt - current time step length
!!  stime - current simulation time
!!
!!***
subroutine Simulation_adjustEvolution(blkcnt, blklst, nstep, dt, stime)

#include "Flash.h"
#include "constants.h"

  use Simulation_data, ONLY : sim_vortSlopeLimit, VORTSL_NONE, VORTSL_VANLEER, &
       VORTSL_MINMOD, VORTSL_MC

  use Grid_interface, ONLY: Grid_getSingleCellVol, &
       Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_fillGuardCells, &
       Grid_releaseBlkPtr, Grid_getSingleCellVol

  implicit none

  integer, intent(in) :: blkcnt
  integer, intent(in) :: blklst(blkcnt)
  integer, intent(in) :: nstep

  real, intent(in) :: dt
  real, intent(in) :: stime

  real :: a
  real :: b
  real :: dx
  real :: dy
  real :: gradx
  real :: grady
  real :: velx
  real :: vely
  real :: delta(MDIM)

  real, pointer :: blkPtr(:,:,:,:)

  integer :: lb
  integer :: i, j, k
  integer :: blkLimitsGC(LOW:HIGH,MDIM)
  integer :: blkLimits(LOW:HIGH,MDIM)

  call Grid_fillGuardCells(CENTER,ALLDIR)

  ! *********************************
  ! *                               *
  ! *     COMPUTE THE VORTICITY     *
  ! *                               *
  ! *********************************

  do lb = 1, blkcnt
     call Grid_getBlkIndexLimits(blklst(lb),blkLimits,blkLimitsGC)
     call Grid_getBlkPtr(blklst(lb), blkPtr)
     call Grid_getDeltas(blklst(lb), delta)

     dx = delta(IAXIS)
     dy = delta(JAXIS)

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              ! In 2D x-y, the voriticy is givien by:
              !   w = khat * (d vy / d x - d vx / dy)

              ! *** COMPUTE VORTICITY WITH SLOPE LIMITERS ***

              ! Compute d vy / dx
              a = blkPtr(VELY_VAR,i+1,j,k)-blkPtr(VELY_VAR,i,j,k)
              b = blkPtr(VELY_VAR,i,j,k)-blkPtr(VELY_VAR,i-1,j,k)

              select case (sim_vortSlopeLimit)
              case (VORTSL_NONE)
                 gradx = nullsl(a,b)/dx                 
              case (VORTSL_VANLEER)
                 gradx = vanLeer(a,b)/dx
              case (VORTSL_MINMOD)
                 gradx = minmod(a,b)/dx
              case (VORTSL_MC)
                 gradx = mc(a,b)/dx
              end select

              ! Compute d vx / dy
              a = blkPtr(VELX_VAR,i,j+1,k)-blkPtr(VELX_VAR,i,j,k)
              b = blkPtr(VELX_VAR,i,j,k)-blkPtr(VELX_VAR,i,j-1,k)

              select case (sim_vortSlopeLimit)
              case (VORTSL_NONE)
                 grady = nullsl(a,b)/dy                
              case (VORTSL_VANLEER)
                 grady = vanLeer(a,b)/dy
              case (VORTSL_MINMOD)
                 grady = minmod(a,b)/dy
              case (VORTSL_MC)
                 grady = mc(a,b)/dy
              end select

              blkPtr(GRDX_VAR,i,j,k) = gradx
              blkPtr(GRDY_VAR,i,j,k) = grady

              blkPtr(VORT_VAR,i,j,k) = gradx - grady
           enddo
        end do
     end do

     call Grid_releaseBlkPtr(blklst(lb), blkPtr)
  end do

contains

  function nullsl(a,b)
    implicit none
    real :: a,b,nullsl
    nullsl = 0.5*(a-b)
  end function nullsl
  
  function minmod(a,b)
    implicit none
    real :: a,b,minmod
    minmod=.5 * (sign(1.,a) + sign(1.,b))*min(abs(a),abs(b))
  end function minmod

  function vanLeer(a,b)
    implicit none
    real :: a,b,vanLeer
    if (a*b <=0.) then
       vanLeer=0.
    else
       vanLeer=2.*a*b/(a+b)
    endif
  end function vanLeer

  function mc(a,b)
    implicit none
    real :: a,b,mc
    mc = (sign(1.,a)+sign(1.,b))*min(abs(a),.25*abs(a+b),abs(b))
  end function mc

  
end subroutine Simulation_adjustEvolution
