!Called from Grid_init
!Allocate and populate the data structure that holds all information
!about the single solid body.

#include "constants.h"
#include "Flash.h"

Subroutine gr_sbInit()
  use gr_sbData, ONLY : gr_sbBodyInfo, gr_sbNumBodies, &
       gr_sbPtNumX, gr_sbPtNumY, gr_sbPtNumZ, gr_sbDebug
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none

  gr_sbNumBodies = 1
  allocate(gr_sbBodyInfo(gr_sbNumBodies))

  gr_sbBodyInfo(1) % boundBox(:,:) = 0.0
  call RuntimeParameters_get("sb_xmin", gr_sbBodyInfo(1) % boundBox(LOW,IAXIS))
  call RuntimeParameters_get("sb_xmax", gr_sbBodyInfo(1) % boundBox(HIGH,IAXIS))
  call RuntimeParameters_get("sb_ptNumX", gr_sbPtNumX)

  if (NDIM >= 2) then
     call RuntimeParameters_get("sb_ymin", gr_sbBodyInfo(1) % boundBox(LOW,JAXIS))
     call RuntimeParameters_get("sb_ymax", gr_sbBodyInfo(1) % boundBox(HIGH,JAXIS))
     call RuntimeParameters_get("sb_ptNumY", gr_sbPtNumY)
  else
     gr_sbPtNumY = 1
  end if

  if (NDIM == 3) then
     call RuntimeParameters_get("sb_zmin", gr_sbBodyInfo(1) % boundBox(LOW,KAXIS))
     call RuntimeParameters_get("sb_zmax", gr_sbBodyInfo(1) % boundBox(HIGH,KAXIS))
     call RuntimeParameters_get("sb_ptNumZ", gr_sbPtNumZ)
  else
     gr_sbPtNumZ = 1
  end if

  call RuntimeParameters_get("sb_debug", gr_sbDebug)
End Subroutine gr_sbInit
