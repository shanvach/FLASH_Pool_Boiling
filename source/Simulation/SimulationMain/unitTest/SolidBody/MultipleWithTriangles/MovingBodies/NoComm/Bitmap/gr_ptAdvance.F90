!!****if* source/Simulation/SimulationMain/unitTest/SolidBody/MultipleWithTriangles/MovingBodies/NoComm/Bitmap/gr_ptAdvance
!!
!! NAME
!!
!!  gr_ptAdvance
!!
!! SYNOPSIS
!!
!!  call gr_ptAdvance(real(in)   :: dtOld,
!!                         real(in)   :: dtNew)
!!
!! DESCRIPTION
!!
!!  Time advancement routine for the body.
!!  
!! ARGUMENTS
!!
!!   dtOld -- not used in this first-order scheme
!!   dtNew -- current time increment

#include "constants.h"
#include "Flash.h"

subroutine gr_ptAdvance (dtOld,dtNew)

  use gr_sbData, ONLY : gr_sbBodyInfo, NumVertices, gr_sbNumBodies
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none

  real, INTENT(in)  :: dtOld, dtNew
  integer       :: i, b
  real          :: jumpx,jumpy,jumpz
  real :: pd_xmin, pd_xmax, pd_ymin, pd_ymax, pd_zmin, pd_zmax

  call RuntimeParameters_get("xmin", pd_xmin)
  call RuntimeParameters_get("xmax", pd_xmax)
  jumpx = dtNew
  if (NDIM >= 2) then
     call RuntimeParameters_get("ymin", pd_ymin)
     call RuntimeParameters_get("ymax", pd_ymax)
     jumpy = dtNew
  else
     jumpy = 0.0
  end if
  if (NDIM == 3) then
     call RuntimeParameters_get("zmin", pd_zmin)
     call RuntimeParameters_get("zmax", pd_zmax)
     jumpz = dtNew
  else
     jumpz = 0.0
  endif
  ! Update the particle positions.
  do b = 1, gr_sbNumBodies

     if (gr_sbBodyInfo(b) % myPE == gr_sbBodyInfo(b) % bodyMaster) then

     NumVertices = gr_sbBodyInfo(b) % NumVertices

     do i = 1, NumVertices
        gr_sbBodyInfo(b) % xb(i) = gr_sbBodyInfo(b) % xb(i) + jumpx
        if ((gr_sbBodyInfo(b) % xb(i) < pd_xmin) .or. (gr_sbBodyInfo(b) % xb(i) > pd_xmax)) then
!           print *, "advnace:x particle outside domain",gr_sbBodyInfo(b) % xb(i)
           gr_sbBodyInfo(b) % xb(i) = NONEXISTENT
        endif
        if (NDIM > 1) then
           gr_sbBodyInfo(b) % yb(i) = gr_sbBodyInfo(b) % yb(i) + jumpy
           if ((gr_sbBodyInfo(b) % yb(i) < pd_ymin) .or. (gr_sbBodyInfo(b) % yb(i) > pd_ymax)) then
!              print *, gr_sbBodyInfo(b) % yb(i)
              gr_sbBodyInfo(b) % yb(i) = NONEXISTENT
           endif
        endif
        if (NDIM > 2) then
           gr_sbBodyInfo(b) % zb(i) = gr_sbBodyInfo(b) % zb(i) + jumpz
           if ((gr_sbBodyInfo(b) % zb(i) < pd_zmin) .or. (gr_sbBodyInfo(b) % zb(i) > pd_zmax)) then
!              print *, gr_sbBodyInfo(b) % zb(i)
              gr_sbBodyInfo(b) % zb(i) = NONEXISTENT
           endif
        endif
     enddo

     endif

  enddo
  return
end subroutine gr_ptAdvance


