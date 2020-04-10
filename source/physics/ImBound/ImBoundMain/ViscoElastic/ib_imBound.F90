!!source/physics/ImBound/ImBoundMain/ViscoElastic
!!
!!
!! NAME
!!
!!  ib_imBound
!!
!!
!! SYNOPSIS
!!
!!  ib_imBound (integer(IN) :: blockCount,
!!             integer(IN) :: blockList(blockCount)
!!             real(IN)    :: timeEndAdv
!!             real(IN)    :: dt)
!!
!!
!! DESCRIPTION
!!
!! ARGUMENTS
!!
!!  blockCount - the number of blocks in blockList
!!  blockList  - array holding local IDs of blocks on which to advance
!!  timeEndAdv - time level at the end of step
!!  dt         - timestep
!!
!!***

subroutine ib_imBound( blockCount, blockList, timeEndAdv, dt)

#include "Flash.h"
#include "ImBound.h"
#include "constants.h"
#include "IncompNS.h"

  ! Modules Use:
#ifdef FLASH_GRID_PARAMESH
  use physicaldata, ONLY : force_consistency,        &
                           interp_mask_unk_res,      &
                           interp_mask_facex_res,    &
                           interp_mask_facey_res,    &
                           interp_mask_facez_res,    &
                           interp_mask_unk,      &
                           interp_mask_facex,    &
                           interp_mask_facey,    &
                           interp_mask_facez
  use workspace, ONLY :    interp_mask_work                           
#endif    

  use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
                             Grid_getListOfBlocks, &
                             Grid_getDeltas,         &
                             Grid_getBlkBC,          &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkIndexLimits, &
                             Grid_fillGuardCells,    &
                             Grid_putFluxData,       &
                             Grid_getFluxData,       &
                             Grid_conserveFluxes,    &
                             Grid_conserveField,     &
                             Grid_updateRefinement,  &
                             Grid_solvePoisson, Grid_getBlkBoundBox, Grid_getBlkCenterCoords

  use gr_interface, ONLY : gr_findMean

  use Grid_Data, ONLY : gr_domainBC 
    
  use Timers_interface, ONLY : Timers_start, Timers_stop

  use Driver_interface, only : Driver_getNStep

  use ib_viscoElastic_interface, only: ib_levelset_linearprojection, ib_levelset_constantprojection, &
                                       ib_dynamic_grid_directional_derivative, ib_solid_stress, ib_ustar_solid, & 
                                       ib_redistance_PM, ib_dynamic_grid_retain_inside, ib_dynamic_grid_normal_vector, & 
                                       ib_solid_interface_advection

  use ib_viscoElastic_data

  use IncompNS_data, only: ins_meshMe

  implicit none

  include "Flash_mpi.h"


  !! ---- Argument List ----------------------------------
  integer, INTENT(INOUT) :: blockCount
  integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList !blockCount
  real,    INTENT(IN) :: timeEndAdv,dt
  !! -----------------------------------------------------

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  real, dimension(2,MDIM) :: boundBox


  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)
            
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

  integer :: lb,blockID,ii,jj,kk,ierr,i,j,k

  integer :: sx,sy,sz,ex,ey,ez

  real dtdxdz,dtdydz,dtdxdy
  
  integer TA(2),count_rate
  real*8  ET

  integer TAIB(2),count_rateIB
  real*8  ETIB

  real maxfp,minfp,maxflb,minflb

  real bsize(MDIM),coord(MDIM)
  integer datasize(MDIM)
  !! define x,y,z coordinate
  !real xcell(MDIM),ycell(MDIM),zcell(MDIM)

  integer nxc, nyc, nzc
  real del(MDIM)

  logical :: gridChanged 

  character(len=6) :: IndNStep
  integer :: NStep
  integer :: step


  if(ins_meshME .eq. MASTER_PE) print *,"Applying IB visco elastic forcing"

!------1: Loop through multiple blocks on a processor
!--------------------calculate components of solid strain tensor---------

  do lb = 1,blockCount
     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC) 

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,facezData,FACEZ)

     call ib_solid_stress(solnData(LMDA_VAR,:,:,:),&
                          solnData(LMDX_VAR,:,:,:),&
                          solnData(LMDY_VAR,:,:,:),&
                          solnData(LMS1_VAR,:,:,:),&
                          solnData(LMS2_VAR,:,:,:),&
                          solnData(LMS3_VAR,:,:,:),&
                          solnData(LMS4_VAR,:,:,:),&
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                          blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
                          del(DIR_X),del(DIR_Y),del(DIR_Z))

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  enddo

  !! Guard Cell Mask
  !gcMask = .FALSE.

  !! BC fill for cell center variables
  !gcMask(LMS1_VAR) = .TRUE.
  !gcMask(LMS2_VAR) = .TRUE.
  !gcMask(LMS3_VAR) = .TRUE.
  !gcMask(LMS4_VAR) = .TRUE.

  !! Fill guard cells
  !call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
  !     maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)

  !------2: Loop through multiple blocks on a processor
  !--------------------update ustar with stress term added---------

  do lb = 1,blockCount
     blockID = blockList(lb)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

     ! Get Blocks internal limits indexes:
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC) 

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,facezData,FACEZ)

     call ib_ustar_solid(facexData(VELC_FACE_VAR,:,:,:),&
                         faceyData(VELC_FACE_VAR,:,:,:),&
                         solnData(XMUS_VAR,:,:,:),&
                         solnData(LMS1_VAR,:,:,:),&
                         solnData(LMS2_VAR,:,:,:),&
                         solnData(LMS3_VAR,:,:,:),&
                         solnData(LMS4_VAR,:,:,:),&
                         blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                         blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                         blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
                         del(DIR_X),del(DIR_Y),del(DIR_Z),dt)


     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  enddo


  ! Guard Cell Mask
  gcMask = .FALSE.

  ! BC fill for face center variables
  gcMask(NUNK_VARS+VELC_FACE_VAR) = .TRUE.
  gcMask(NUNK_VARS+1*NFACE_VARS+VELC_FACE_VAR) = .TRUE.
#if NDIM == 3
  gcMask(NUNK_VARS+2*NFACE_VARS+VELC_FACE_VAR) = .TRUE.
#endif

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           
 
end subroutine ib_imBound
