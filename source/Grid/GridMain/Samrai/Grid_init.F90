!!****if* source/Grid/GridMain/Samrai/Grid_init
!!
!! NAME
!!  Grid_init
!!
!! SYNOPSIS
!!  #include "Flash.h"
!!
!!  Grid_init(integer(IN) :: myPE,
!!            integer(IN) :: numProcs,
!!            restart(IN) :: restart)
!!
!! DESCRIPTION
!!  Initialize the Grid
!!
!!***

subroutine Grid_init(myPE, numProcs, restart)

  use Grid_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Simulation_interface, ONLY : Simulation_mapIntToStr
implicit none
#include "Flash.h"
#include "constants.h"

  integer, intent(IN) :: myPE,numProcs
  logical, intent(IN) :: restart
  integer :: i

  character(len=4),dimension(UNK_VARS_BEGIN:UNK_VARS_END) :: gr_unkLabels

  call RuntimeParameters_get('xmin', gr_physDomainMinMax(1))
  call RuntimeParameters_get('xmax', gr_physDomainMinMax(4))
  call RuntimeParameters_get('ymin', gr_physDomainMinMax(2))
  call RuntimeParameters_get('ymax', gr_physDomainMinMax(5))
  call RuntimeParameters_get('zmin', gr_physDomainMinMax(3))
  call RuntimeParameters_get('zmax' ,gr_physDomainMinMax(6))

  gr_indexDomainMinMax = 1
  call RuntimeParameters_get("iGridSize", gr_indexDomainMinMax(2))
  call RuntimeParameters_get("jGridSize", gr_indexDomainMinMax(4))
  call RuntimeParameters_get("kGridSize", gr_indexDomainMinMax(6))   


  call RuntimeParameters_get("iguard", gr_iguard)
  call RuntimeParameters_get("jguard", gr_jguard)
  call RuntimeParameters_get("kguard", gr_kguard)

  call RuntimeParameters_get("refine_ratio", gr_refineRatio)

  call RuntimeParameters_get("imaxPatchSize",gr_indexMaxPatchSize(1))
  call RuntimeParameters_get("jmaxPatchSize",gr_indexMaxPatchSize(2))
  call RuntimeParameters_get("kmaxPatchSize",gr_indexMaxPatchSize(3))

  call RuntimeParameters_get("iminPatchSize",gr_indexMinPatchSize(1))
  call RuntimeParameters_get("jminPatchSize",gr_indexMinPatchSize(2))
  call RuntimeParameters_get("kminPatchSize",gr_indexMinPatchSize(3))

  call RuntimeParameters_get("efficiencyTolerance", gr_effTolerance)
  call RuntimeParameters_get("combineEfficiency", gr_combineEff)
  call RuntimeParameters_get("lrefine_max", lrefine_max)

  do i = UNK_VARS_BEGIN,UNK_VARS_END
     call Simulation_mapIntToStr(i,gr_unkLabels,MAPBLOCK_UNK)
  end do

  call Samrai_init(gr_unkLabels, NUNK_VARS, guard, gr_physDomainMinMax,&
                   gr_indexDomainMinMax,procGrid, lrefine_max, refine_ratio,&
                   gr_indexMaxPatchSize, gr_indexMinPatchSize, effTolerance,&
                   combineEff, argc, argv )



!! Dev : Anshu this is for completeness, but doesn't
!! mean anything yet
!!$  call RuntimeParameters_get("xl_boundary_type", xl_bcString)
!!$  call RuntimeParameters_get("xr_boundary_type", xr_bcString)
!!$  call RuntimeParameters_get("yl_boundary_type", yl_bcString)
!!$  call RuntimeParameters_get("yr_boundary_type", yr_bcString)
!!$  call RuntimeParameters_get("zl_boundary_type", zl_bcString)
!!$  call RuntimeParameters_get("zr_boundary_type", zr_bcString)
  call RuntimeParameters_get("geometry", str_geometry)

!!$  call RuntimeParameters_get("priority_dir1",priorityOneDir)
!!$  call RuntimeParameters_get("priority_dir2",priorityTwoDir)
!!$
!!$  call RuntimeParameters_mapStrToInt(xl_bcString,il_bcType)
!!$  call RuntimeParameters_mapStrToInt(xr_bcString,ir_bcType)
!!$  call RuntimeParameters_mapStrToInt(yl_bcString,jl_bcType)
!!$  call RuntimeParameters_mapStrToInt(yr_bcString,jr_bcType)
!!$  call RuntimeParameters_mapStrToInt(zl_bcString,kl_bcType)
!!$  call RuntimeParameters_mapStrToInt(zr_bcString,kr_bcType)



end subroutine Grid_init
