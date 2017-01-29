!*******************************************************************************

!  Routine:     mg_restrict()

!  Description: Restrict a multigrid variable from one mesh level to the
!               next coarser level.  Assume that child blocks on the "from"
!               level and leaf blocks on the "to" level contain valid data.
!               On output, all blocks on the "to" level should then contain
!               valid data.

!  Parameters:  level       Level to restrict from.
!               ifrom       Index of the variable to restrict from.
!               ito         Index of the variable to restrict into.


subroutine gr_mgRestrictFaces (level, ifrom, ito)

!===============================================================================

  use mg_common, ONLY: nodetype_save, ili,jli,kli,iui,jui,kui,ile, jle, kle, iue, jue, kue

  use Grid_data, ONLY : gr_meshMe,gr_meshNumProcs

  use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                                Grid_getBlkPtr,       &
                                Grid_releaseBlkPtr,   &
                                Grid_fillGuardCells,  &
                                Grid_getBlkCenterCoords, &
                                Grid_getDeltas,       &
                                Grid_getBlkBoundBox


  use tree, only : nodetype,lrefine

  use workspace , ONLY: work, interp_mask_work

  use physicaldata, ONLY : force_consistency,        &
                           interp_mask_unk_res,      &
                           interp_mask_facex_res,    &
                           interp_mask_facey_res,    &
                           interp_mask_facez_res,    &
                           interp_mask_unk,      &
                           interp_mask_facex,    &
                           interp_mask_facey,    &
                           interp_mask_facez

  !- kpd - Added
  use paramesh_interfaces, ONLY: amr_restrict

  use Driver_data, ONLY: dr_nstep

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: level, ifrom, ito

  integer :: indexPlotNumber, PTNumber, ierr, intval
  integer :: lb, lnblocks, i,j,k 
  real    :: normf, normt
  real, pointer, dimension(:,:,:,:),save :: facevarx , facevary, facevarz, &
                                            facexData, faceyData, facezData

  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  real, dimension(2,MDIM) :: boundBox
  real coord(MDIM),xcell,ycell
  real del(MDIM),bsize(MDIM)


!===============================================================================

!               The PARAMESH library's restriction and boundary-update routines
!               only operate on all arrays at once or on the "work" array.
!               So we copy the input variable to work(), restrict work(),
!               then copy to the output variable.  However, we have to trick
!               the mesh package into thinking all of the "from"-level blocks
!               are leaf nodes, so that the supplied restriction routine can
!               do its work.

!               Temporarily alter node type arrays.

  call Grid_getLocalNumBlks(lnblocks)

!===============================================================================
!===============================================================================
!===============================================================================


  do lb = 1, lnblocks

     !--------------------------------------------
     !- kpd - TEST LOOP
     !--------------------------------------------

          ! Point to blocks face vars:
          call Grid_getBlkPtr(lb,facevarx ,FACEX)
          call Grid_getBlkPtr(lb,facexData,FACEX)
          call Grid_getBlkPtr(lb,facevary ,FACEY)          
          call Grid_getBlkPtr(lb,faceyData,FACEY)

          ! Release pointers:
          call Grid_releaseBlkPtr(lb,facevarx ,FACEX)
          call Grid_releaseBlkPtr(lb,facevary ,FACEY)
          call Grid_releaseBlkPtr(lb,facexData,FACEX)
          call Grid_releaseBlkPtr(lb,faceyData,FACEY)
  enddo

!===============================================================================
!===============================================================================
!===============================================================================
  !- kpd - Set interpolation (restriction) order, affects weighting...
  !        -> 2 = Lagrange (-1/8,6/8,3/8) polynomial
  !        -> 1 = Equal (1/2,1/2) weighting
  intval = 1
  !intval = 2
!!  interp_mask_unk = intval;   interp_mask_unk_res = intval;
!!  interp_mask_work= intval;
  interp_mask_facex = intval; interp_mask_facex_res = intval;
  interp_mask_facey = intval; interp_mask_facey_res = intval;
  interp_mask_facez = intval; interp_mask_facez_res = intval;

  !print*,"Entering amr_restrictFaces"
  !--------------------------------------
  !- kpd - Now perform the restriction...
  !*******************************************************
  call amr_restrictFaces (gr_meshMe, 1, 0, .false., level)
  !call amr_restrictFaces (gr_meshMe, 1, 0, .true., level)
  !*******************************************************
  !print*,"Leaving amr_restrictFaces"

  !- kpd - This line is in the original, doesn't seem to change soln...
  !call amr_get_new_nodetypes (gr_meshNumProcs, gr_meshMe, level-1)

!===============================================================================
!===============================================================================
!===============================================================================

  do lb = 1, lnblocks

     !--------------------------------------------
     !- kpd - If the block is on the coarser level
     !--------------------------------------------

          ! Point to blocks face vars:
          call Grid_getBlkPtr(lb,facevarx ,FACEX)
          call Grid_getBlkPtr(lb,facevary ,FACEY)
          call Grid_getBlkPtr(lb,facexData,FACEX)
          call Grid_getBlkPtr(lb,faceyData,FACEY)

          ! Release pointers:
          call Grid_releaseBlkPtr(lb,facevarx ,FACEX)
          call Grid_releaseBlkPtr(lb,facevary ,FACEY)
          call Grid_releaseBlkPtr(lb,facexData,FACEX)
          call Grid_releaseBlkPtr(lb,faceyData,FACEY)
  enddo

!===============================================================================
!===============================================================================
!===============================================================================


          !------------------------------------------------
          !- kpd - Fill Guard Cells for the mixture density
          !------------------------------------------------
          gcMask                                       = .FALSE.
          gcMask(NUNK_VARS+MGW8_FACE_VAR)              = .TRUE.    ! rhoX
          gcMask(NUNK_VARS+1*NFACE_VARS+MGW8_FACE_VAR) = .TRUE.    ! rhoY
#if NDIM == 3
          gcMask(NUNK_VARS+2*NFACE_VARS+MGW8_FACE_VAR) = .TRUE.    ! rhoZ
#endif

  intval = 1
  !intval = 2
  interp_mask_facex = intval; interp_mask_facex_res = intval;
  interp_mask_facey = intval; interp_mask_facey_res = intval;
  interp_mask_facez = intval; interp_mask_facez_res = intval;

          call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
               maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)

   !print*,"Leaving GC Fill amr_restrictFaces"
!===============================================================================
!===============================================================================
!===============================================================================

  !-----------------
  !- kpd - Test Loop
  !-----------------

  do lb = 1, lnblocks

     !--------------------------------------------
     !- kpd - If the block is on the coarser level
     !--------------------------------------------
     !if (lrefine(lb) == level-1) then

          ! Point to blocks face vars:
          call Grid_getBlkPtr(lb,facevarx ,FACEX)
          call Grid_getBlkPtr(lb,facevary ,FACEY)
          call Grid_getBlkPtr(lb,facexData,FACEX)
          call Grid_getBlkPtr(lb,faceyData,FACEY)
#if NDIM == 3
          call Grid_getBlkPtr(lb,facevarz ,FACEZ)
          call Grid_getBlkPtr(lb,facezData,FACEZ)
#endif

!          !- kpd - Get block limits (including GCs)
!          call Grid_getBlkIndexLimits(lb,blkLimits,blkLimitsGC)
!
          call Grid_getBlkCenterCoords(lb,coord)
          call Grid_getBlkBoundBox(lb,boundBox)
          call Grid_getDeltas(lb,del)

          bsize(:) = boundBox(2,:) - boundBox(1,:)

          do k=1,kue         !1
             do j=1,jue      !blkLimitsGC(HIGH,JAXIS)
                do i=1,iue   !blkLimitsGC(HIGH,IAXIS)

                      xcell = coord(IAXIS) - bsize(IAXIS)/2.0 +   &
                              real(i - NGUARD - 1)*del(IAXIS) +   &
                              0.5*del(IAXIS)

                      ycell  = coord(JAXIS) - bsize(JAXIS)/2.0 +  &
                              real(j - NGUARD - 1)*del(JAXIS)  +  &
                              0.5*del(JAXIS)

!                   if (dr_nstep .eq. 14) then
                      !if (lb .eq. 18) print*,"Final",MGW8_FACE_VAR,ito,"--:>",lb,i,j,facevarx(ito,i,j,k),facevary(ito,i,j,k)
!                     !print*,"Final",gr_meshMe,lb,coord(1),coord(2),i,j,xcell,ycell,facevarx(ito,i,j,k),facevary(ito,i,j,k)
!                   end if

                    !print*,"Final",MGW8_FACE_VAR,ito,i,j,facevarx(ito,i,j,k),facevary(ito,i,j,k)

                end do
             end do
          end do

          ! Release pointers:
          call Grid_releaseBlkPtr(lb,facevarx ,FACEX)
          call Grid_releaseBlkPtr(lb,facevary ,FACEY)
          call Grid_releaseBlkPtr(lb,facexData,FACEX)
          call Grid_releaseBlkPtr(lb,faceyData,FACEY)
#if NDIM == 3
          call Grid_releaseBlkPtr(lb,facevarz ,FACEZ)
          call Grid_releaseBlkPtr(lb,facezData,FACEZ)
#endif

     !endif
  enddo

  !print*,"Leaving gr_mgRestrictFaces All Together"

!===============================================================================

  return
end subroutine gr_mgRestrictFaces 
