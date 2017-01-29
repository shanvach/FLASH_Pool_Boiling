!*******************************************************************************

!  Routine:     mg_gr_mgAssignFden()

!  Description: Initializes a multigrid variable on a particular
!               mesh level.

!  Parameters:  level       Level to zero on.
!               ivar        Variable to initialize.
!               leaf_only   If nonzero, only zero leaf-node blocks on this
!                           level.


subroutine gr_mgAssignFden (level, ivar, leaf_only)

!===============================================================================

use mg_common, ONLY: nodetype_save, ili, jli, kli, iui, jui, kui,&
     ile, jle, kle, iue, jue, kue

use tree, only : lrefine
use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                              Grid_getBlkPtr,       &
                              Grid_releaseBlkPtr,   &
                              Grid_fillGuardCells

  use physicaldata, ONLY : force_consistency,        &
                           interp_mask_unk_res,      &
                           interp_mask_facex_res,    &
                           interp_mask_facey_res,    &
                           interp_mask_facez_res,    &
                           interp_mask_unk,      &
                           interp_mask_facex,    &
                           interp_mask_facey,    &
                           interp_mask_facez

implicit none
#include "constants.h"
#include "Flash.h"

integer, intent(in) :: level, ivar, leaf_only
integer :: lnblocks, blockID, intval
integer :: lb, i, j, k
real, pointer, dimension(:,:,:,:) :: facevarx,facevary,facevarz
real, pointer, dimension(:,:,:,:) :: facexData, faceyData, facezData
real :: Mdens

logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)

!===============================================================================

call Grid_getLocalNumBlks(lnblocks)

if (leaf_only == 0) then

  do lb = 1, lnblocks

    !-------------------------------------
    !- kpd - Loop throught ALL leaf blocks
    !-------------------------------------
!    if (lrefine(lb) == level .OR. &
!        ((lrefine(lb) .lt. level)  .and. (nodetype_save(lb) == 1)) ) then

      ! Point to blocks center vars:
      call Grid_getBlkPtr(lb,facexData,FACEX)
      call Grid_getBlkPtr(lb,faceyData,FACEY)
#if NDIM == 3
      call Grid_getBlkPtr(lb,facezData,FACEZ)
#endif

      !print*,"Assigning density to",lb,lrefine(lb),nodetype_save(lb)

      !- kpd - Include Guard Cells In Loop
      do k=1,kue
       do j=1,jue
          do i=1,iue


             Mdens = facexData(RH1F_FACE_VAR,i,j,k) + &
                     facexData(RH2F_FACE_VAR,i,j,k)

             facexData(ivar,i,j,k) = Mdens
             !!!facexData(ivar,i,j,k) = 1.0

             Mdens = faceyData(RH1F_FACE_VAR,i,j,k) + &
                     faceyData(RH2F_FACE_VAR,i,j,k)

             faceyData(ivar,i,j,k) = Mdens
             !!!faceyData(ivar,i,j,k) = 1.0

#if NDIM == 3
             Mdens = facezData(RH1F_FACE_VAR,i,j,k) + &
                     facezData(RH2F_FACE_VAR,i,j,k)

             facezData(ivar,i,j,k) = Mdens
#endif

!             if (lb .eq. 41 .OR. lb .eq. 29 .OR. lb .eq. 40) print*,"AssignFden",ivar,lb,i,j,facexData(ivar,i,j,k),faceyData(ivar,i,j,k)

             !if (lrefine(lb) == 1) print*,"Level1 Density",lb,i,j,facexData(ivar,i,j,k),faceyData(ivar,i,j,k)
             !print*,"gr_mgAssign Density",lb,i,j,facexData(ivar,i,j,k),faceyData(ivar,i,j,k)

          enddo
        enddo
      enddo

      ! Release pointers:
      call Grid_releaseBlkPtr(lb,facexData,FACEX)
      call Grid_releaseBlkPtr(lb,faceyData,FACEY)
#if NDIM == 3
      call Grid_releaseBlkPtr(lb,facezData,FACEZ)
#endif

!    endif
  enddo

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

!- kpd - Shouldn't need to do this, b/c RH1F_FACE_VAR & RH2F_FACE_VAR should 
!           be appropriately filled already!

!          !------------------------------------------------
!          !- kpd - Fill Guard Cells for the mixture density
!          !         Note: MGW8_FACE_VAR = ivar
!          !------------------------------------------------
!          gcMask                                       = .FALSE.
!          gcMask(NUNK_VARS+MGW8_FACE_VAR)              = .TRUE.    ! rhoX
!          gcMask(NUNK_VARS+1*NFACE_VARS+MGW8_FACE_VAR) = .TRUE.    ! rhoY
!#if NDIM == 3
!          gcMask(NUNK_VARS+2*NFACE_VARS+MGW8_FACE_VAR) = .TRUE.    ! rhoZ
!#endif
!
!          !-------------------------------------------------------------------
!          !- kpd - Set interpolation (restriction) order, affects weighting...
!          !        -> 2 = Lagrange (-1/8,6/8,3/8) polynomial - UNBOUNDED!
!          !        -> 1 = Equal (1/2,1/2) weighting - BOUNDED!
!          !-------------------------------------------------------------------
!          intval = 1
!          interp_mask_unk = intval;   interp_mask_unk_res = intval;
!          interp_mask_facex = intval; interp_mask_facex_res = intval;
!          interp_mask_facey = intval; interp_mask_facey_res = intval;
!          interp_mask_facez = intval; interp_mask_facez_res = intval;
!
!          call Grid_fillGuardCells(CENTER_FACES,ALLDIR,&
!               maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

else

   print*,"Error in AssignFdens, the leaf_only flag is wrong!"

endif

!===============================================================================

return
end
