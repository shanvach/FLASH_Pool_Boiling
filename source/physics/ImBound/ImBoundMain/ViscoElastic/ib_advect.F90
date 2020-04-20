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

subroutine ib_advect( blockCount, blockList, timeEndAdv, dt)

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
                                       ib_lsRedistance, ib_advectWENO3, ib_solid_interface_advection, ib_fluid_props

  use ib_viscoElastic_data

  use Grid_data, only: gr_meshMe

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

  real :: minCellDiag, lsT, lsDT

  integer :: max_lsit, maxiter, intval

  integer :: time_projection(2)
  real*8  :: elapsed_time

  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: lmda_old


  CALL SYSTEM_CLOCK(time_projection(1),count_rate)

  max_lsit = 3
  maxiter = 200

#define SOLID_RECONSTRUCTION 0
#define BC_OPTIMIZE 1

#ifdef FLASH_GRID_PARAMESH
    !intval = 1
    intval = 2
    interp_mask_unk = intval;   interp_mask_unk_res = intval;
    interp_mask_work= intval;
#endif

  if(gr_meshMe .eq. MASTER_PE) print *,"Entering IB level set advection" 

  !------1: Loop through multiple blocks on a processor
  !--------------------WENO3 advection-----------

  ! WENO3 advection for LMDX, LMDY
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

     call ib_advectWENO3(solnData(LMDA_VAR,:,:,:), &
                         solnData(LMDX_VAR,:,:,:), &
                         facexData(VELC_FACE_VAR,:,:,:), &
                         faceyData(VELC_FACE_VAR,:,:,:), &
                         dt, &
                         del(DIR_X), &
                         del(DIR_Y), &
                         blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
#if BC_OPTIMIZE == 1
                         blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),.TRUE.)
#else
                         blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),.FALSE.)
#endif

     call ib_advectWENO3(solnData(LMDA_VAR,:,:,:), &
                         solnData(LMDY_VAR,:,:,:), &
                         facexData(VELC_FACE_VAR,:,:,:), &
                         faceyData(VELC_FACE_VAR,:,:,:), &
                         dt, &
                         del(DIR_X), &
                         del(DIR_Y), &
                         blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
#if BC_OPTIMIZE == 1
                         blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),.TRUE.)
#else
                         blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),.FALSE.)
#endif

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  enddo

  ! Guard Cell Mask
  gcMask = .FALSE.
  
  ! BC fill for cell center variables
  gcMask(LMDX_VAR) = .TRUE.
  gcMask(LMDY_VAR) = .TRUE.
  
  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           

  !------2: Loop through multiple blocks on a processor
  !--------------------calculate normal vectors---------

  !!!calculate normal vector
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

     call ib_dynamic_grid_normal_vector(solnData(LMDA_VAR,:,:,:),&
                                  solnData(NMLX_VAR,:,:,:),&
                                  solnData(NMLY_VAR,:,:,:),&
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

  ! Guard Cell Mask
  gcMask = .FALSE.
  
  ! BC fill for cell center variables
  gcMask(NMLX_VAR) = .TRUE.
  gcMask(NMLY_VAR) = .TRUE.
  
  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           
 
  !------3: Loop through multiple blocks on a processor
  !--------------------dynamic grid projection for X grid---------

  !!!define directional derivative
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

     call ib_dynamic_grid_directional_derivative(solnData(LMDA_VAR,:,:,:),&
                                  solnData(LMDX_VAR,:,:,:),&
                                  solnData(NMLX_VAR,:,:,:),&
                                  solnData(NMLY_VAR,:,:,:),&
                                  solnData(DDSN_VAR,:,:,:),&
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

  ! Guard Cell Mask
  gcMask = .FALSE.
  
  ! BC fill for cell center variables
  gcMask(DDSN_VAR) = .TRUE.
   
  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           

  do step = 1, maxiter !projection step can be changed
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

     call ib_levelset_constantprojection(solnData(LMDA_VAR,:,:,:),&
                                         solnData(DDSN_VAR,:,:,:),&
                                         solnData(NMLX_VAR,:,:,:),&
                                         solnData(NMLY_VAR,:,:,:),&
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

#if BC_OPTIMIZE == 1
  enddo !end of iteration

  ! Guard Cell Mask
  gcMask = .FALSE.

  ! BC fill for cell center variables
  gcMask(DDSN_VAR) = .TRUE.

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           
#else
  ! Guard Cell Mask
  gcMask = .FALSE.

  ! BC fill for cell center variables
  gcMask(DDSN_VAR) = .TRUE.

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           

  enddo !end of iteration
#endif

  do step = 1, maxiter !projection step can be changed
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

     call ib_levelset_linearprojection(solnData(LMDA_VAR,:,:,:),&
                                       solnData(LMDX_VAR,:,:,:),&
                                       solnData(DDSN_VAR,:,:,:),&
                                       solnData(NMLX_VAR,:,:,:),&
                                       solnData(NMLY_VAR,:,:,:),&
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

#if BC_OPTIMIZE == 1
  enddo !end of iteration

  ! Guard Cell Mask
  gcMask = .FALSE.

  ! BC fill for cell center variables
  gcMask(LMDX_VAR) = .TRUE.

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           
#else
  ! Guard Cell Mask
  gcMask = .FALSE.

  ! BC fill for cell center variables
  gcMask(LMDX_VAR) = .TRUE.

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           

  enddo !end of iteration
#endif
  !!!end of linear extrapolation of X grid

  !------4: Loop through multiple blocks on a processor
  !--------------------dynamic grid projection for Y grid---------
 
  !!!define directional derivative
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

     call ib_dynamic_grid_directional_derivative(solnData(LMDA_VAR,:,:,:),&
                                  solnData(LMDY_VAR,:,:,:),&
                                  solnData(NMLX_VAR,:,:,:),&
                                  solnData(NMLY_VAR,:,:,:),&
                                  solnData(DDSN_VAR,:,:,:),&
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

  ! Guard Cell Mask
  gcMask = .FALSE.
  
  ! BC fill for cell center variables
  gcMask(DDSN_VAR) = .TRUE.
  
  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           

  do step = 1, maxiter !projection step can be changed
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

     call ib_levelset_constantprojection(solnData(LMDA_VAR,:,:,:),&
                                         solnData(DDSN_VAR,:,:,:),&
                                         solnData(NMLX_VAR,:,:,:),&
                                         solnData(NMLY_VAR,:,:,:),&
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

#if BC_OPTIMIZE == 1
  enddo !end of iteration

  ! Guard Cell Mask
  gcMask = .FALSE.

  ! BC fill for cell center variables
  gcMask(DDSN_VAR) = .TRUE.

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           
#else
  ! Guard Cell Mask
  gcMask = .FALSE.

  ! BC fill for cell center variables
  gcMask(DDSN_VAR) = .TRUE.

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           

  enddo !end of iteration
#endif

  do step = 1, maxiter !projection step can be changed
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

     call ib_levelset_linearprojection(solnData(LMDA_VAR,:,:,:),&
                                  solnData(LMDY_VAR,:,:,:),&
                                  solnData(DDSN_VAR,:,:,:),&
                                  solnData(NMLX_VAR,:,:,:),&
                                  solnData(NMLY_VAR,:,:,:),&
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

#if BC_OPTIMIZE == 1
  enddo !end of iteration

  ! Guard Cell Mask
  gcMask = .FALSE.

  ! BC fill for cell center variables
  gcMask(LMDY_VAR) = .TRUE.

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           
#else
  ! Guard Cell Mask
  gcMask = .FALSE.

  ! BC fill for cell center variables
  gcMask(LMDY_VAR) = .TRUE.

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           

  enddo !end of iteration
#endif

  !------5: Loop through multiple blocks on a processor
  !--------------------calculate the level set of solid interface either using advection or reconstruction with X-Y grid---------

#if SOLID_RECONSTRUCTION == 1

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

     !ib_solid_interface_advection needs to be changed for different initial interface profile
     call ib_solid_interface_advection(solnData(LMDA_VAR,:,:,:),&
                                       solnData(LMDX_VAR,:,:,:),&
                                       solnData(LMDY_VAR,:,:,:),&
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

  ! Guard Cell Mask
  gcMask = .FALSE.
   
  ! BC fill for cell center variables
  gcMask(LMDA_VAR) = .TRUE.
  
  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           

#else

  ! WENO3 advection for LMDA
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

     call ib_advectWENO3(solnData(LMDA_VAR,:,:,:), &
                         solnData(LMDA_VAR,:,:,:), &
                         facexData(VELC_FACE_VAR,:,:,:), &
                         faceyData(VELC_FACE_VAR,:,:,:), &
                         dt, &
                         del(DIR_X), &
                         del(DIR_Y), &
                         blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                         blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),.TRUE.)

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  enddo

  ! Guard Cell Mask
  gcMask = .FALSE.
  
  ! BC fill for cell center variables
  gcMask(LMDA_VAR) = .TRUE.
 
  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           


  ! Re-initialization for LMDA
  do ii=1,max_lsit
 
  lsT = 0.0

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

     minCellDiag = SQRT(del(DIR_X)**2.+del(DIR_Y)**2.)
     lsDT = minCellDiag/2.0d0
     if (ii .eq. max_lsit .AND. lb .eq. 1 .AND. gr_meshMe .eq. 0) then
       print*,"Level Set Initialization Iteration # ",ii,minCellDiag,lsDT
     end if

     if (ii.eq.1) lmda_old = solnData(LMDA_VAR,:,:,:)

     call ib_lsRedistance(solnData(LMDA_VAR,:,:,:), &
                          del(DIR_X),del(DIR_Y),  &
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS), &
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS), &
                          lmda_old, lsDT, blockID,minCellDiag)

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  enddo

  ! Guard Cell Mask
  gcMask = .FALSE.
  
  ! BC fill for cell center variables
  gcMask(LMDA_VAR) = .TRUE.
 
  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           

  lsT = lsT + lsDT
  enddo

#endif

  !------6: Loop through multiple blocks on a processor
  !--------------------assign fluid properties---------

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

     call ib_fluid_props(solnData(LMDA_VAR,:,:,:),&
                         ib_vis_rho1, ib_vis_rho2,&
                         ib_vis_xmu1, ib_vis_xmu2,& 
                         ib_vis_xmus,             &
                         solnData(XRHO_VAR,:,:,:),&
                         solnData(XMUF_VAR,:,:,:),&
                         solnData(XMUS_VAR,:,:,:),&
                         blockID,                 &
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

  ! Guard Cell Mask
  gcMask = .FALSE.
  
  ! BC fill for cell center variables
  gcMask(XRHO_VAR) = .TRUE.
  gcMask(XMUF_VAR) = .TRUE.
  gcMask(XMUS_VAR) = .TRUE.

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           
 
  CALL SYSTEM_CLOCK(time_projection(2),count_rate)
  elapsed_time=REAL(time_projection(2)-time_projection(1),8)/count_rate
  if (gr_meshMe .eq. MASTER_PE)  write(*,*) 'Total Level Set Advection Time =', elapsed_time

end subroutine ib_advect
