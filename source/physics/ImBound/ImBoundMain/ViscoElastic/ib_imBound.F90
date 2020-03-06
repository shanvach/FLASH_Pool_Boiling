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

  use ib_viscoElastic_interface

  use ib_viscoElastic_data

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

  integer nxc, nyc, nzc
  real del(MDIM)

  logical :: gridChanged 

  character(len=6) :: IndNStep
  integer :: NStep
  integer :: step

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

     ! Place function calls here - Below is sample function call
     !call ib_sample_function_call(facexData(VARF_FACE_VAR,:,:,:),&
     !                             faceyData(VARF_FACE_VAR,:,:,:),&
     !                             facezData(VARF_FACE_VAR,:,:,:),&
     !                              solnData(VARC_VAR,:,:,:),&
     !                  blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
     !                  blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
     !                  blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS),&
     !                  del(DIR_X),del(DIR_Y),del(DIR_Z))

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


  ! Guard Cell Mask
  gcMask = .FALSE.

  ! BC fill for cell center variables
  !gcMask(VARC_VAR) = .TRUE.
  gcMask(LMDX_VAR) = .TRUE.
  gcMask(LMDY_VAR) = .TRUE.
  gcMask(LMS1_VAR) = .TRUE.
  gcMask(LMS2_VAR) = .TRUE.
  gcMask(LMS3_VAR) = .TRUE.
  gcMask(LMS4_VAR) = .TRUE.

!  ! BC fill for face center variables
!  gcMask(NUNK_VARS+VARF_FACE_VAR) = .TRUE.
!  gcMask(NUNK_VARS+1*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#if NDIM == 3
!  gcMask(NUNK_VARS+2*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#endif

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)

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
  !gcMask(VARC_VAR) = .TRUE.
  gcMask(MUSF_VAR) = .TRUE.
  gcMask(LMS1_VAR) = .TRUE.
  gcMask(LMS2_VAR) = .TRUE.
  gcMask(LMS3_VAR) = .TRUE.
  gcMask(LMS4_VAR) = .TRUE.

  ! BC fill for face center variables
  gcMask(NUNK_VARS+VELC_FACE_VAR) = .TRUE.
  gcMask(NUNK_VARS+1*NFACE_VARS+VELC_FACE_VAR) = .TRUE.
#if NDIM == 3
  gcMask(NUNK_VARS+2*NFACE_VARS+VELC_FACE_VAR) = .TRUE.
#endif

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           
 
  !------3 should be called after LMDX&LMDY advection------
  !------3: Loop through multiple blocks on a processor
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


     !ib_solid_interface_advection(sd,sX,sY,&
     !                                ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz)

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
  !gcMask(VARC_VAR) = .TRUE.
  gcMask(LMDA_VAR) = .TRUE.
  gcMask(LMDX_VAR) = .TRUE.
  gcMask(LMDY_VAR) = .TRUE.

!  ! BC fill for face center variables
!  gcMask(NUNK_VARS+VARF_FACE_VAR) = .TRUE.
!  gcMask(NUNK_VARS+1*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#if NDIM == 3
!  gcMask(NUNK_VARS+2*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#endif

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           
 
  !------4: Loop through multiple blocks on a processor
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


     !subroutine ib_redistance_PM(s,ib_vis_rho1, ib_vis_rho2, ib_vis_xmu1, ib_vis_xmu2,&
     !                                   ib_vis_xmus, rho, xmu, xms,blockID,&
     !                                   ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz)
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


  ! Guard Cell Mask
  gcMask = .FALSE.

  ! BC fill for cell center variables
  !gcMask(VARC_VAR) = .TRUE.
  gcMask(LMDA_VAR) = .TRUE.
  gcMask(XRHO_VAR) = .TRUE.
  gcMask(XMUF_VAR) = .TRUE.
  gcMask(XMUS_VAR) = .TRUE.

!  ! BC fill for face center variables
!  gcMask(NUNK_VARS+VARF_FACE_VAR) = .TRUE.
!  gcMask(NUNK_VARS+1*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#if NDIM == 3
!  gcMask(NUNK_VARS+2*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#endif

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           
 
  !------5: Loop through multiple blocks on a processor
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
  !gcMask(VARC_VAR) = .TRUE.
  gcMask(LMDA_VAR) = .TRUE.
  gcMask(LMDX_VAR) = .TRUE.
  gcMask(ADFX_VAR) = .TRUE.
  gcMask(ADFY_VAR) = .TRUE.

!  ! BC fill for face center variables
!  gcMask(NUNK_VARS+VARF_FACE_VAR) = .TRUE.
!  gcMask(NUNK_VARS+1*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#if NDIM == 3
!  gcMask(NUNK_VARS+2*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#endif

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
  !gcMask(VARC_VAR) = .TRUE.
  gcMask(LMDA_VAR) = .TRUE.
  gcMask(LMDX_VAR) = .TRUE.
  gcMask(ADFX_VAR) = .TRUE.
  gcMask(ADFY_VAR) = .TRUE.
  gcMask(DDSN_VAR) = .TRUE.

!  ! BC fill for face center variables
!  gcMask(NUNK_VARS+VARF_FACE_VAR) = .TRUE.
!  gcMask(NUNK_VARS+1*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#if NDIM == 3
!  gcMask(NUNK_VARS+2*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#endif

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           
 
  !!!constant projection of directional derivative
  solnData(D0SN_VAR,:,:,:) = solnData(DDSN_VAR,:,:,:)
  do step = 1, 200 !projection step can be changed
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



     call ib_levelset_constantprojection(solnData(D0SN_VAR,:,:,:),&
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
  !gcMask(VARC_VAR) = .TRUE.
  gcMask(ADFX_VAR) = .TRUE.
  gcMask(ADFY_VAR) = .TRUE.
  gcMask(D0SN_VAR) = .TRUE.

!  ! BC fill for face center variables
!  gcMask(NUNK_VARS+VARF_FACE_VAR) = .TRUE.
!  gcMask(NUNK_VARS+1*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#if NDIM == 3
!  gcMask(NUNK_VARS+2*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#endif

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           
 
! retain level set inside solid. Update level set outside solid 
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



     call ib_dynamic_grid_retain_inside(solnData(LMDA_VAR,:,:,:),&
                                  solnData(D0SN_VAR,:,:,:),&
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
  !gcMask(VARC_VAR) = .TRUE.
  gcMask(LMDA_VAR) = .TRUE.
  gcMask(D0SN_VAR) = .TRUE.
  gcMask(DDSN_VAR) = .TRUE.

!  ! BC fill for face center variables
!  gcMask(NUNK_VARS+VARF_FACE_VAR) = .TRUE.
!  gcMask(NUNK_VARS+1*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#if NDIM == 3
!  gcMask(NUNK_VARS+2*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#endif

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           
 
! retain level set inside solid. Update level set outside solid

enddo
solnData(DDSN_VAR,:,:,:) = solnData(D0SN_VAR,:,:,:)
!!!end of constant projection of directional derivative

!!!linear extrapolation of X grid
        solnData(LM0X_VAR,:,:,:) = solnData(LMDX_VAR,:,:,:)
        do step = 1, 200 !projection step can be changed
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



     call ib_levelset_linearprojection(solnData(LM0X_VAR,:,:,:),&
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


  ! Guard Cell Mask
  gcMask = .FALSE.

  ! BC fill for cell center variables
  !gcMask(VARC_VAR) = .TRUE.
  gcMask(ADFX_VAR) = .TRUE.
  gcMask(ADFY_VAR) = .TRUE.
  gcMask(LM0X_VAR) = .TRUE.
  gcMask(DDSN_VAR) = .TRUE.

!  ! BC fill for face center variables
!  gcMask(NUNK_VARS+VARF_FACE_VAR) = .TRUE.
!  gcMask(NUNK_VARS+1*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#if NDIM == 3
!  gcMask(NUNK_VARS+2*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#endif

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           

! ! retain level set inside solid. Update level set outside solid 
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



     call ib_dynamic_grid_retain_inside(solnData(LMDA_VAR,:,:,:),&
                                  solnData(LM0X_VAR,:,:,:),&
                                  solnData(LMDX_VAR,:,:,:),&
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
  !gcMask(VARC_VAR) = .TRUE.
  gcMask(LMDA_VAR) = .TRUE.
  gcMask(LM0X_VAR) = .TRUE.
  gcMask(LMDX_VAR) = .TRUE.

!  ! BC fill for face center variables
!  gcMask(NUNK_VARS+VARF_FACE_VAR) = .TRUE.
!  gcMask(NUNK_VARS+1*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#if NDIM == 3
!  gcMask(NUNK_VARS+2*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#endif

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           
 
! ! end of retain level set inside solid. Update level set outside solid 

        enddo
        solnData(LMDX_VAR,:,:,:) = solnData(LM0X_VAR,:,:,:)
!!!end of linear extrapolation of X grid

  !------6: Loop through multiple blocks on a processor
  !--------------------dynamic grid projection for Y grid---------
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
  !gcMask(VARC_VAR) = .TRUE.
  gcMask(LMDA_VAR) = .TRUE.
  gcMask(LMDY_VAR) = .TRUE.
  gcMask(ADFX_VAR) = .TRUE.
  gcMask(ADFY_VAR) = .TRUE.

!  ! BC fill for face center variables
!  gcMask(NUNK_VARS+VARF_FACE_VAR) = .TRUE.
!  gcMask(NUNK_VARS+1*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#if NDIM == 3
!  gcMask(NUNK_VARS+2*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#endif

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
  !gcMask(VARC_VAR) = .TRUE.
  gcMask(LMDA_VAR) = .TRUE.
  gcMask(LMDY_VAR) = .TRUE.
  gcMask(ADFX_VAR) = .TRUE.
  gcMask(ADFY_VAR) = .TRUE.
  gcMask(DDSN_VAR) = .TRUE.

!  ! BC fill for face center variables
!  gcMask(NUNK_VARS+VARF_FACE_VAR) = .TRUE.
!  gcMask(NUNK_VARS+1*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#if NDIM == 3
!  gcMask(NUNK_VARS+2*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#endif

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           
 
  !!!constant projection of directional derivative
  solnData(D0SN_VAR,:,:,:) = solnData(DDSN_VAR,:,:,:)
  do step = 1, 200 !projection step can be changed
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



     call ib_levelset_constantprojection(solnData(D0SN_VAR,:,:,:),&
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
  !gcMask(VARC_VAR) = .TRUE.
  gcMask(ADFX_VAR) = .TRUE.
  gcMask(ADFY_VAR) = .TRUE.
  gcMask(D0SN_VAR) = .TRUE.

!  ! BC fill for face center variables
!  gcMask(NUNK_VARS+VARF_FACE_VAR) = .TRUE.
!  gcMask(NUNK_VARS+1*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#if NDIM == 3
!  gcMask(NUNK_VARS+2*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#endif

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           
 
! retain level set inside solid. Update level set outside solid 
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



     call ib_dynamic_grid_retain_inside(solnData(LMDA_VAR,:,:,:),&
                                  solnData(D0SN_VAR,:,:,:),&
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
  !gcMask(VARC_VAR) = .TRUE.
  gcMask(LMDA_VAR) = .TRUE.
  gcMask(D0SN_VAR) = .TRUE.
  gcMask(DDSN_VAR) = .TRUE.

!  ! BC fill for face center variables
!  gcMask(NUNK_VARS+VARF_FACE_VAR) = .TRUE.
!  gcMask(NUNK_VARS+1*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#if NDIM == 3
!  gcMask(NUNK_VARS+2*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#endif

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           
 
! retain level set inside solid. Update level set outside solid

enddo
solnData(DDSN_VAR,:,:,:) = solnData(D0SN_VAR,:,:,:)
!!!end of constant projection of directional derivative

!!!linear extrapolation of Y grid
        solnData(LM0Y_VAR,:,:,:) = solnData(LMDY_VAR,:,:,:)
        do step = 1, 200 !projection step can be changed
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



     call ib_levelset_linearprojection(solnData(LM0Y_VAR,:,:,:),&
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


  ! Guard Cell Mask
  gcMask = .FALSE.

  ! BC fill for cell center variables
  !gcMask(VARC_VAR) = .TRUE.
  gcMask(ADFX_VAR) = .TRUE.
  gcMask(ADFY_VAR) = .TRUE.
  gcMask(LM0Y_VAR) = .TRUE.
  gcMask(DDSN_VAR) = .TRUE.

!  ! BC fill for face center variables
!  gcMask(NUNK_VARS+VARF_FACE_VAR) = .TRUE.
!  gcMask(NUNK_VARS+1*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#if NDIM == 3
!  gcMask(NUNK_VARS+2*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#endif

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           

! ! retain level set inside solid. Update level set outside solid 
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



     call ib_dynamic_grid_retain_inside(solnData(LMDA_VAR,:,:,:),&
                                  solnData(LM0Y_VAR,:,:,:),&
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
  !gcMask(VARC_VAR) = .TRUE.
  gcMask(LMDA_VAR) = .TRUE.
  gcMask(LM0Y_VAR) = .TRUE.
  gcMask(LMDY_VAR) = .TRUE.

!  ! BC fill for face center variables
!  gcMask(NUNK_VARS+VARF_FACE_VAR) = .TRUE.
!  gcMask(NUNK_VARS+1*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#if NDIM == 3
!  gcMask(NUNK_VARS+2*NFACE_VARS+VARF_FACE_VAR) = .TRUE.
!#endif

  ! Fill guard cells
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)           
 
! ! end of retain level set inside solid. Update level set outside solid 

        enddo
        solnData(LMDY_VAR,:,:,:) = solnData(LM0Y_VAR,:,:,:)
!!!end of linear extrapolation of X grid

 


end subroutine ib_imBound
