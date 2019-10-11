subroutine ib_advect (blockCount,blockList,dt)

#include "Flash.h"
#include "ImBound.h"
#include "constants.h"

#define IB_REDISTANCE

  ! Modules Use:
#ifdef FLASH_GRID_PARAMESH
     use physicaldata, ONLY : interp_mask_unk_res,      &
                              interp_mask_facex_res,    &
                              interp_mask_facey_res,    &
                              interp_mask_facez_res,    &
                              interp_mask_unk,      &
                              interp_mask_facex,    &
                              interp_mask_facey,    &
                             interp_mask_facez
     use workspace, ONLY :    interp_mask_work
#endif    

     use Grid_interface, ONLY : Grid_getListOfBlocks, &
                                Grid_getDeltas,         &
                                Grid_getBlkBC,          &
                                Grid_getBlkPtr,         &
                                Grid_releaseBlkPtr,     &
                                Grid_getBlkIndexLimits, &
                                Grid_fillGuardCells,    &
                                Grid_getBlkBoundBox,Grid_getBlkCenterCoords


      use ib_lset_interface, only: ib_advectWENO3, ib_advectWENO3_3D, ib_lsRedistance, &
                                   ib_lsRedistance_3D

      use Driver_data, ONLY : dr_nstep, dr_simTime

      use IncompNS_data, ONLY : ins_meshMe,ins_alfa,ins_gravX,ins_gravY,ins_invRe,ins_gravZ

      implicit none

#include "IncompNS.h"

  include "Flash_mpi.h"

      integer, INTENT(INOUT) :: blockCount
      integer, INTENT(INOUT), dimension(MAXBLOCKS) :: blockList
      real,    INTENT(IN) :: dt

      integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
      real, dimension(2,MDIM) :: boundBox
      real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData
      integer :: lb,blockID,ii,jj,kk,ierr,i,j,k
      real bsize(MDIM),coord(MDIM)
      real del(MDIM),xcell,ycell, ycellp,zcell,rc,xcellp,zcellp

      integer :: listofBlocks(MAXBLOCKS)
      integer :: count
      integer :: intval

      logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)
      real :: lsDT,lsT,minCellDiag
  
      integer :: max_lsit

      real :: amp, omega, pi, offset, phase

      amp = 0.5
      omega = 0.1
      pi = acos(-1.0)
      offset = 0.5
      phase = pi

      max_lsit = 3

      do lb = 1,blockCount

        blockID = blockList(lb)

        call Grid_getBlkBoundBox(blockId,boundBox)

        bsize(:) = boundBox(2,:) - boundBox(1,:)

        call Grid_getBlkCenterCoords(blockId,coord)

        ! Get blocks dx, dy ,dz:
        call Grid_getDeltas(blockID,del)

        ! Get Blocks internal limits indexes:
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

        ! Point to blocks center and face vars:
        call Grid_getBlkPtr(blockID,solnData,CENTER)
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)
        call Grid_getBlkPtr(blockID,facezData,FACEZ)

        facexData(VELB_FACE_VAR,:,:,:) =   0.0

        !faceyData(VELB_FACE_VAR,:,:,:) =  0.0
        faceyData(VELB_FACE_VAR,:,:,:) = -1.0
        !faceyData(VELB_FACE_VAR,:,:,:) =  2*pi*omega*amp*sin(2*pi*omega*dr_simTime);

#if NDIM == 2
        call ib_advectWENO3(solnData(LMDA_VAR,:,:,:), &
                          facexData(VELB_FACE_VAR,:,:,:), &
                          faceyData(VELB_FACE_VAR,:,:,:), &
                          dt, &
                          del(DIR_X), &
                          del(DIR_Y), &
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS))
#endif

#if NDIM == 3
        facezData(VELB_FACE_VAR,:,:,:) = 0.0
        call ib_advectWENO3_3D(solnData(LMDA_VAR,:,:,:), &
                          facexData(VELB_FACE_VAR,:,:,:), &
                          faceyData(VELB_FACE_VAR,:,:,:), &
                          facezData(VELB_FACE_VAR,:,:,:), &
                          dt, &
                          del(DIR_X), &
                          del(DIR_Y), &
                          del(DIR_Z), &
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS),&
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS),&
                          blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS))
#endif

     ! Release pointers:
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

    end do

!--------------------------------------------------------------------
    gcMask = .FALSE.
    gcMask(LMDA_VAR) = .TRUE.

    call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask,selectBlockType=ACTIVE_BLKS)
!--------------------------------------------------------------------

#ifdef IB_REDISTANCE
   do ii = 1,max_lsit

     !------------------------------
     !- kpd - Level set redistancing 
     !------------------------------

     !lsDT = dt
     lsT  = 0.0

     do lb = 1,blockCount
        blockID = blockList(lb)

         call Grid_getBlkBoundBox(blockId,boundBox)
         bsize(:) = boundBox(2,:) - boundBox(1,:)

        call Grid_getBlkCenterCoords(blockId,coord)

        ! Get blocks dx, dy ,dz:
        call Grid_getDeltas(blockID,del)

        ! Get Blocks internal limits indexes:
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

        ! Point to blocks center and face vars:
        call Grid_getBlkPtr(blockID,solnData,CENTER)
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)
        call Grid_getBlkPtr(blockID,facezData,FACEZ)

#if NDIM == 3
        lsDT = MIN(10.0*dt,0.001)
        minCellDiag =SQRT((SQRT(del(DIR_X)**2.+del(DIR_Y)**2.))**2.+del(DIR_Z)**2.)
        if ( ii .eq. max_lsit .AND. lb .eq. 1 .AND. ins_meshMe .eq. 0) then
           print*,"IB Level Set Initialization Iteration # ",ii,minCellDiag,lsDT
        end if

        if (ii.eq.1) solnData(AAJUNK_VAR,:,:,:) = solnData(LMDA_VAR,:,:,:)

        call ib_lsRedistance_3D(solnData(LMDA_VAR,:,:,:), &
                          facexData(VELB_FACE_VAR,:,:,:), &
                          faceyData(VELB_FACE_VAR,:,:,:), &
                          facezData(VELB_FACE_VAR,:,:,:), &
                          del(DIR_X),del(DIR_Y),del(DIR_Z),  &
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS), &
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS), &
                          blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS), &
                          solnData(AAJUNK_VAR,:,:,:), lsDT, minCellDiag )

#elif NDIM == 2

        minCellDiag = SQRT(del(DIR_X)**2.+del(DIR_Y)**2.)
        lsDT = minCellDiag/2.0d0
        if ( ii .eq. max_lsit .AND. lb .eq. 1 .AND. ins_meshMe .eq. 0) then
           print*,"IB Level Set Initialization Iteration # ",ii,minCellDiag,lsDT
        end if

        if (ii.eq.1) solnData(AAJUNK_VAR,:,:,:) = solnData(LMDA_VAR,:,:,:)

        call ib_lsRedistance(solnData(LMDA_VAR,:,:,:), &
                          facexData(VELB_FACE_VAR,:,:,:),  &
                          faceyData(VELB_FACE_VAR,:,:,:),  &
                          del(DIR_X),del(DIR_Y),  &
                          blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS), &
                          blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS), &
                          solnData(AAJUNK_VAR,:,:,:), lsDT, blockID,minCellDiag)

#endif
        ! Release pointers:
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
     enddo

    gcMask = .FALSE.
    gcMask(LMDA_VAR) = .TRUE.
    call Grid_fillGuardCells(CENTER,ALLDIR,&
       maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)

    lsT = lsT + lsDT

   end do
#endif

#if NDIM == 2
  do lb=1,blockCount

        ! Loop through all the blocks
        blockID = blockList(lb)

        ! Get the block Id and boundaries
        call Grid_getBlkBoundBox(blockId,boundBox)

        ! Size of the block (far side - near side)
        bsize(:) = boundBox(2,:) - boundBox(1,:)

        ! Get co-ordinates of the block:
        call Grid_getBlkCenterCoords(blockId,coord)

        ! Get the delta x and y for the block:
        call Grid_getDeltas(blockID,del)

        ! Get Blocks internal limits indexes:
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

        ! Point to blocks center and face vars:
        call Grid_getBlkPtr(blockID,solnData,CENTER)
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)
        call Grid_getBlkPtr(blockID,facezData,FACEZ)

        k = 1
        do j=2,blkLimitsGC(HIGH,JAXIS)-1
            do i=2,blkLimitsGC(HIGH,IAXIS)-1

               solnData(NMLX_VAR,i,j,k) = -((solnData(LMDA_VAR,i+1,j,k) - solnData(LMDA_VAR,i-1,j,k))/2*del(IAXIS))/&
                                      sqrt(((solnData(LMDA_VAR,i+1,j,k) - solnData(LMDA_VAR,i-1,j,k))/2*del(IAXIS))**2+&
                                           ((solnData(LMDA_VAR,i,j+1,k) - solnData(LMDA_VAR,i,j-1,k))/2*del(JAXIS))**2)

               solnData(NMLY_VAR,i,j,k) = -((solnData(LMDA_VAR,i,j+1,k) - solnData(LMDA_VAR,i,j-1,k))/2*del(IAXIS))/&
                                      sqrt(((solnData(LMDA_VAR,i+1,j,k) - solnData(LMDA_VAR,i-1,j,k))/2*del(IAXIS))**2+&
                                           ((solnData(LMDA_VAR,i,j+1,k) - solnData(LMDA_VAR,i,j-1,k))/2*del(JAXIS))**2)


            end do
        end do

        do j=2,blkLimitsGC(HIGH,JAXIS)-1
            do i=2,blkLimitsGC(HIGH,IAXIS)-1
               solnData(TNGY_VAR,i,j,k) = ((solnData(LMDA_VAR,i+1,j,k) - solnData(LMDA_VAR,i-1,j,k))/2*del(IAXIS))/&
                                     sqrt(((solnData(LMDA_VAR,i+1,j,k) - solnData(LMDA_VAR,i-1,j,k))/2*del(IAXIS))**2+&
                                          ((solnData(LMDA_VAR,i,j+1,k) - solnData(LMDA_VAR,i,j-1,k))/2*del(JAXIS))**2)

               solnData(TNGX_VAR,i,j,k) = -((solnData(LMDA_VAR,i,j+1,k) - solnData(LMDA_VAR,i,j-1,k))/2*del(IAXIS))/&
                                      sqrt(((solnData(LMDA_VAR,i+1,j,k) - solnData(LMDA_VAR,i-1,j,k))/2*del(IAXIS))**2+&
                                           ((solnData(LMDA_VAR,i,j+1,k) - solnData(LMDA_VAR,i,j-1,k))/2*del(JAXIS))**2)

            end do
        end do

        ! Release pointers:
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  end do
#endif

#if NDIM == 3
  do lb=1,blockCount

        ! Loop through all the blocks
        blockID = blockList(lb)

        ! Get the block Id and boundaries
        call Grid_getBlkBoundBox(blockId,boundBox)

        ! Size of the block (far side - near side)
        bsize(:) = boundBox(2,:) - boundBox(1,:)

        ! Get co-ordinates of the block:
        call Grid_getBlkCenterCoords(blockId,coord)

        ! Get the delta x and y for the block:
        call Grid_getDeltas(blockID,del)

        ! Get Blocks internal limits indexes:
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

        ! Point to blocks center and face vars:
        call Grid_getBlkPtr(blockID,solnData,CENTER)
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)
        call Grid_getBlkPtr(blockID,facezData,FACEZ)

        do k=2,blkLimitsGC(HIGH,KAXIS)-1
        do j=2,blkLimitsGC(HIGH,JAXIS)-1
        do i=2,blkLimitsGC(HIGH,IAXIS)-1

           solnData(NMLX_VAR,i,j,k) = -((solnData(LMDA_VAR,i+1,j,k) - solnData(LMDA_VAR,i-1,j,k))/2*del(IAXIS))/&
                                      sqrt(((solnData(LMDA_VAR,i+1,j,k) - solnData(LMDA_VAR,i-1,j,k))/2*del(IAXIS))**2+&
                                           ((solnData(LMDA_VAR,i,j+1,k) - solnData(LMDA_VAR,i,j-1,k))/2*del(JAXIS))**2+&
                                           ((solnData(LMDA_VAR,i,j,k+1) - solnData(LMDA_VAR,i,j,k-1))/2*del(KAXIS))**2)

           solnData(NMLY_VAR,i,j,k) = -((solnData(LMDA_VAR,i,j+1,k) - solnData(LMDA_VAR,i,j-1,k))/2*del(JAXIS))/&
                                      sqrt(((solnData(LMDA_VAR,i+1,j,k) - solnData(LMDA_VAR,i-1,j,k))/2*del(IAXIS))**2+&
                                           ((solnData(LMDA_VAR,i,j+1,k) - solnData(LMDA_VAR,i,j-1,k))/2*del(JAXIS))**2+&
                                           ((solnData(LMDA_VAR,i,j,k+1) - solnData(LMDA_VAR,i,j,k-1))/2*del(KAXIS))**2)

           solnData(NMLZ_VAR,i,j,k) = -((solnData(LMDA_VAR,i,j,k+1) - solnData(LMDA_VAR,i,j,k-1))/2*del(KAXIS))/&
                                      sqrt(((solnData(LMDA_VAR,i+1,j,k) - solnData(LMDA_VAR,i-1,j,k))/2*del(IAXIS))**2+&
                                           ((solnData(LMDA_VAR,i,j+1,k) - solnData(LMDA_VAR,i,j-1,k))/2*del(JAXIS))**2+&
                                           ((solnData(LMDA_VAR,i,j,k+1) - solnData(LMDA_VAR,i,j,k-1))/2*del(KAXIS))**2)

        end do
        end do
        end do

        ! Release pointers:
        call Grid_releaseBlkPtr(blockID,solnData,CENTER)
        call Grid_releaseBlkPtr(blockID,facexData,FACEX)
        call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
        call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  end do
#endif

end subroutine ib_advect
