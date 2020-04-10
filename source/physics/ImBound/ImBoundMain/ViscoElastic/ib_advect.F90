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
                                       ib_lsRedistance, ib_advectWENO3, ib_solid_interface_advection

  use ib_viscoElastic_data

  use IncompNS_data, only: ins_meshME
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

  integer :: max_lsit, maxiter

  integer :: time_projection(2)
  real*8  :: elapsed_time

  max_lsit = 3
  maxiter = 60

  CALL SYSTEM_CLOCK(time_projection(1),count_rate)

  if(ins_meshME .eq. MASTER_PE) print *,"Entering IB level set advection" 

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

     call ib_advectWENO3(solnData(LMDX_VAR,:,:,:), &
                         facexData(VELC_FACE_VAR,:,:,:), &
                         faceyData(VELC_FACE_VAR,:,:,:), &
                         dt, &
                         del(DIR_X), &
                         del(DIR_Y), &
                         blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                         blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS))

     call ib_advectWENO3(solnData(LMDY_VAR,:,:,:), &
                         facexData(VELC_FACE_VAR,:,:,:), &
                         faceyData(VELC_FACE_VAR,:,:,:), &
                         dt, &
                         del(DIR_X), &
                         del(DIR_Y), &
                         blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                         blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS))

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

  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           

  !------2: Loop through multiple blocks on a processor
  !--------------------dynamic grid projection for X grid---------

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
                                  solnData(LMDX_VAR,:,:,:),&
                                  solnData(ADFX_VAR,:,:,:),&
                                  solnData(ADFY_VAR,:,:,:),&
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
  gcMask(ADFX_VAR) = .TRUE.
  gcMask(ADFY_VAR) = .TRUE.
  
  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           
 
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
                                  solnData(ADFX_VAR,:,:,:),&
                                  solnData(ADFY_VAR,:,:,:),&
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
                                         solnData(ADFX_VAR,:,:,:),&
                                         solnData(ADFY_VAR,:,:,:),&
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
  enddo !end of iteration

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

     call ib_levelset_linearprojection(solnData(LMDA_VAR,:,:,:),&
                                       solnData(LMDX_VAR,:,:,:),&
                                       solnData(DDSN_VAR,:,:,:),&
                                       solnData(ADFX_VAR,:,:,:),&
                                       solnData(ADFY_VAR,:,:,:),&
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
  enddo !end of iteration

  ! Guard Cell Mask
  gcMask = .FALSE.

  ! BC fill for cell center variables
  gcMask(LMDX_VAR) = .TRUE.

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           

  !!!end of linear extrapolation of X grid

!  !------3: Loop through multiple blocks on a processor
!  !--------------------dynamic grid projection for Y grid---------

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
                                  solnData(LMDY_VAR,:,:,:),&
                                  solnData(ADFX_VAR,:,:,:),&
                                  solnData(ADFY_VAR,:,:,:),&
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
  gcMask(ADFX_VAR) = .TRUE.
  gcMask(ADFY_VAR) = .TRUE.

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           
 
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
                                  solnData(ADFX_VAR,:,:,:),&
                                  solnData(ADFY_VAR,:,:,:),&
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
                                         solnData(ADFX_VAR,:,:,:),&
                                         solnData(ADFY_VAR,:,:,:),&
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
  enddo !end of iteration

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

     call ib_levelset_linearprojection(solnData(LMDA_VAR,:,:,:),&
                                  solnData(LMDY_VAR,:,:,:),&
                                  solnData(DDSN_VAR,:,:,:),&
                                  solnData(ADFX_VAR,:,:,:),&
                                  solnData(ADFY_VAR,:,:,:),&
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
  enddo !end of iteration

  ! Guard Cell Mask
  gcMask = .FALSE.

  ! BC fill for cell center variables
  gcMask(LMDY_VAR) = .TRUE.

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           

  !------4: Loop through multiple blocks on a processor
  !--------------------calculate the level set of interface using X and Y grid---------

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
 
  !------5: Loop through multiple blocks on a processor
  !--------------------redistancing the level set of interface using projection method---------

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

     call ib_redistance_PM(solnData(LMDA_VAR,:,:,:),&
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

  !! Guard Cell Mask
  !gcMask = .FALSE.

  !! BC fill for cell center variables
  !gcMask(XRHO_VAR) = .TRUE.
  !gcMask(XMUF_VAR) = .TRUE.
  !gcMask(XMUS_VAR) = .TRUE.
  !!gcMask(LMDA_VAR) = .TRUE.

  !! Fill guard cells
  !call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
  !     maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           
 
  CALL SYSTEM_CLOCK(time_projection(2),count_rate)
  elapsed_time=REAL(time_projection(2)-time_projection(1),8)/count_rate
  if (ins_meshMe .eq. MASTER_PE)  write(*,*) 'Total Level Set Advection Time =', elapsed_time

end subroutine ib_advect
