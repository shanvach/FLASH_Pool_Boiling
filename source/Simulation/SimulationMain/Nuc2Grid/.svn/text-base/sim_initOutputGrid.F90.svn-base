!!****if* source/Simulation/SimulationMain/Nuc2Grid/sim_initOutputGrid
!!
!! NAME
!!
!!  sim_initOutputGrid
!!
!!
!! SYNOPSIS
!!
!!  call sim_initOutputGrid()
!!
!!
!! DESCRIPTION
!!
!!  Initializes all the parameters needed for the output grid of the
!!  nucOutput-to-radTrans converter.
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!  sim_smlRho : 
!!  sim_gam : 
!!***

subroutine sim_initOutputGrid()
  
  use sim_outputGridData

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, RuntimeParameters_mapStrToInt
  use Logfile_interface, ONLY : Logfile_stamp
  use Grid_data, ONLY : gr_nBlockX, gr_nBlockY, gr_nBlockZ  ! added by Min

  implicit none
#include "Flash.h"


  real :: CdXmin, CdXMax, CdYmin, CdYMax, CdZmin, CdZmax
  real :: OgRatioI, OgRatioJ, OgRatioK
  integer :: OgLowerI, OgLowerJ, OgLowerK
  real :: checkGrid

  call Logfile_stamp("initializing for NucOToRT output grid",  &
       "[sim_initOutputGrid]")

  !call RuntimeParameters_get("sim_densityThreshold",sim_densityThreshold)
  call RuntimeParameters_get("radTranDataFile",sim_ogRadTranDataFileName)

!  Determine size of the output domain

!  call RuntimeParameters_get("radTranOutputXMin",xOgMin)
!  call RuntimeParameters_get("radTranOutputXMax",xOgMax)
!  call RuntimeParameters_get("radTranOutputYMin",yOgMin)
!  call RuntimeParameters_get("radTranOutputYMax",yOgMax)
!  call RuntimeParameters_get("radTranOutputZMin",zOgMin)
!  call RuntimeParameters_get("radTranOutputZMax",zOgMax)

  call RuntimeParameters_get("xmin",CdXMin)
  call RuntimeParameters_get("xmax",CdXMax)
  call RuntimeParameters_get("ymin",CdYMin)
  call RuntimeParameters_get("ymax",CdYMax)
  call RuntimeParameters_get("zmin",CdZMin)
  call RuntimeParameters_get("zmax",CdZMax)

  call RuntimeParameters_get("OgRatioI",OgRatioI)
  call RuntimeParameters_get("OgRatioJ",OgRatioJ)
  call RuntimeParameters_get("OgRatioK",OgRatioK)

  xOgMin = CdXmin*OgRatioI
  xOgMax = CdXmax*OgRatioI
  yOgMin = CdYmin*OgRatioJ
  yOgMax = CdYmax*OgRatioJ
  zOgMin = CdZmin*OgRatioK
  zOgMax = CdZmax*OgRatioK

! Control the resolution of the output domain/grid

!  call RuntimeParameters_get("radTranGridSizeI",nIOg)
!  call RuntimeParameters_get("radTranGridSizeJ",nJOg)
!  call RuntimeParameters_get("radTranGridSizeK",nKOg)

  call RuntimeParameters_get("OgLowerI",OgLowerI)
  call RuntimeParameters_get("OgLowerJ",OgLowerJ)
  call RuntimeParameters_get("OgLowerK",OgLowerK)

  call RuntimeParameters_get("radTranOutputNdim",ndimOg)
  call RuntimeParameters_get("radTranOutputGeometry",geometryOgStr)
  call RuntimeParameters_mapStrToInt(geometryOgStr, geometryOg)

  checkGrid = NXB*gr_nBlockX*OgRatioI/OgLowerI
  nIOg = checkGrid
  if (nIOg-checkGrid/=0) then
     print*, 'Wrong, Og I does not coincide in boundary of cells'
     return
  end if

  if(NDIM<2) then
    nJOg = 1
  else 
    checkGrid = NYB*gr_nBlockY*OgRatioJ/OgLowerJ
    nJOg = checkGrid
    if (nJOg-checkGrid/=0) then
        print*, 'Wrong, Og J does not coincide in boundary of cells'
        return
     end if
  end if

  if(NDIM<3) then
    nKOg = 1
  else
    checkGrid = NZB*gr_nBlockZ*OgRatioK/OgLowerK
    nKOg = checkGrid
    if (nKOg-checkGrid/=0) then
        print*, 'Wrong, Og K does not coincide in boundary of cells'
        return
     end if
  end if 

  iOgB  = 1; iOgE  = nIOg
  iOgB0 = 0; iOgE0 = nIOg+1
  if (ndimOg>1) then
     kOg2D = 1
     jOgB  = 1; jOgE  = nJOg
     jOgB0 = 0; jOgE0 = nJOg+1
  else
     kOg2D = 0
     jOgB  = 1; jOgE  = 1
     jOgB0 = 1; jOgE0 = 1
  end if
  if (ndimOg>2) then
     kOg3D = 1
     kOgB  = 1; kOgE  = nKOg
     kOgB0 = 0; kOgE0 = nKOg+1
  else
     kOg3D = 0
     kOgB  = 1; kOgE  = 1
     kOgB0 = 1; kOgE0 = 1
  end if

end subroutine sim_initOutputGrid
